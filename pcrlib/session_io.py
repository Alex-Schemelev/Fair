import os
import shutil
import time

import pandas as pd


def load_dimers_matrix_from_run(run_dir):
    """Подгрузить матрицу димеров с диска (не храним в session.json)."""
    path = os.path.join(run_dir, "primer_dimers_matrix.csv")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    try:
        df = pd.read_csv(path, sep=None, engine="python")
        return df.to_dict(orient="records")
    except Exception:
        return []


def is_path_within_directory(directory, target):
    abs_directory = os.path.abspath(directory)
    abs_target = os.path.abspath(target)
    try:
        return os.path.commonpath([abs_directory, abs_target]) == abs_directory
    except ValueError:
        return False


def safe_extract_zip(zf, dest_dir):
    """Распаковка zip с защитой от Zip Slip."""
    dest_dir = os.path.abspath(dest_dir)
    os.makedirs(dest_dir, exist_ok=True)
    for info in zf.infolist():
        name = info.filename.replace("\\", "/")
        if not name or name.endswith("/"):
            # каталог — создадим при извлечении файлов
            continue
        # отбрасываем абсолютные и parent-пути
        if os.path.isabs(name) or name.startswith("/") or name.startswith("../") or "/../" in f"/{name}/":
            raise ValueError(f"Небезопасный путь в архиве: {info.filename}")
        target = os.path.abspath(os.path.join(dest_dir, name))
        if not is_path_within_directory(dest_dir, target):
            raise ValueError(f"Zip Slip заблокирован: {info.filename}")
    zf.extractall(dest_dir)
    # повторная проверка на случай странностей extractall
    for root, _dirs, files in os.walk(dest_dir):
        for fn in files:
            full = os.path.join(root, fn)
            if not is_path_within_directory(dest_dir, full):
                raise ValueError(f"После распаковки найден файл вне каталога: {full}")


def slim_session_payload(payload):
    """Убрать тяжёлые поля из session.json (они остаются файлами в run_dir)."""
    out = dict(payload or {})
    out.pop("dimers_matrix", None)
    out["dimers_matrix_file"] = "primer_dimers_matrix.csv"
    return out


def cleanup_old_runs(runs_dir, max_age_days=30, keep_latest=20):
    """
    Удалить каталоги прогонов старше max_age_days (по mtime),
    но всегда сохранить keep_latest самых свежих по mtime.
    """
    if not runs_dir or not os.path.isdir(runs_dir):
        return
    try:
        entries = []
        for name in os.listdir(runs_dir):
            path = os.path.join(runs_dir, name)
            if os.path.isdir(path):
                try:
                    mtime = os.path.getmtime(path)
                except Exception:
                    continue
                entries.append((mtime, path))
        if not entries:
            return
        entries.sort(key=lambda x: x[0], reverse=True)
        keep = set(path for _mt, path in entries[: max(0, int(keep_latest))])
        cutoff = time.time() - (float(max_age_days) * 86400.0)
        for mtime, path in entries:
            if path in keep:
                continue
            if mtime < cutoff:
                try:
                    shutil.rmtree(path, ignore_errors=True)
                except Exception:
                    pass
    except Exception:
        pass
