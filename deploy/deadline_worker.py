"""Sync worker that refuses to wait forever on a half-open client.

Gunicorn's ``timeout`` also caps a long PCR run, so it cannot be short.
The previous gthread setup never applied that timeout to a request thread:
a client that connected and stopped sending held ``recv()`` for days, the
accept loop kept sockets until the open-file limit, and the port stopped
answering.
"""

import errno
import ssl
import time

from gunicorn import http
from gunicorn import sock
from gunicorn import util
from gunicorn.workers.sync import SyncWorker

# Silence longer than this between socket reads aborts the request.
IDLE_READ_TIMEOUT = 20
# A trickle of bytes cannot hold a worker longer than this.
TOTAL_READ_TIMEOUT = 120
# A client that does not read the response cannot hold a worker either.
SEND_TIMEOUT = 120


class ReadDeadlineSocket:
    """Proxy that bounds request reads without limiting the PCR computation.

    The computation itself does not touch the socket. Only ``recv`` is on a
    deadline; ``send`` gets its own timeout so a stalled download is dropped.
    """

    def __init__(self, raw, idle_timeout, total_timeout):
        self._sock = raw
        self._idle = idle_timeout
        self._total = total_timeout
        self._started = None

    def recv(self, bufsize, flags=0):
        now = time.monotonic()
        if self._started is None:
            self._started = now
        remaining = self._total - (now - self._started)
        if remaining <= 0:
            raise TimeoutError("request read deadline exceeded")
        self._sock.settimeout(min(self._idle, remaining))
        return self._sock.recv(bufsize, flags)

    def send(self, data, flags=0):
        self._sock.settimeout(SEND_TIMEOUT)
        return self._sock.send(data, flags)

    def sendall(self, data, flags=0):
        self._sock.settimeout(SEND_TIMEOUT)
        return self._sock.sendall(data, flags)

    def close(self):
        self._sock.close()

    def __getattr__(self, name):
        return getattr(self._sock, name)


class DeadlineWorker(SyncWorker):
    def handle(self, listener, client, addr):
        client = ReadDeadlineSocket(client, IDLE_READ_TIMEOUT, TOTAL_READ_TIMEOUT)
        req = None
        try:
            if self.cfg.is_ssl:
                client = sock.ssl_wrap_socket(client, self.cfg)
            parser = http.RequestParser(self.cfg, client, addr)
            req = next(parser)
            self.handle_request(listener, req, client, addr)
        except TimeoutError:
            self.log.warning("Timed out reading request from %s", addr)
        except http.errors.NoMoreData as e:
            self.log.debug("Ignored premature client disconnection. %s", e)
        except StopIteration as e:
            self.log.debug("Closing connection. %s", e)
        except ssl.SSLError as e:
            if e.args[0] == ssl.SSL_ERROR_EOF:
                self.log.debug("ssl connection closed")
                client.close()
            else:
                self.log.debug("Error processing SSL request.")
                self.handle_error(req, client, addr, e)
        except OSError as e:
            if e.errno not in (errno.EPIPE, errno.ECONNRESET, errno.ENOTCONN):
                self.log.exception("Socket error processing request.")
            else:
                self.log.debug("Ignored socket error %s", e.errno)
        except BaseException as e:
            self.handle_error(req, client, addr, e)
        finally:
            util.close(client)
