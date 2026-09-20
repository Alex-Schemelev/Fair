@echo off
setlocal EnableExtensions
cd /d "%~dp0"

set "VENV_PY=%~dp0.venv\Scripts\python.exe"
set "PORTABLE_PY=%~dp0portable_python\python.exe"

if exist "%VENV_PY%" (
  set "PY=%VENV_PY%"
) else if exist "%PORTABLE_PY%" (
  echo [WARN] .venv not found. Using portable_python directly.
  echo Run setup_env.bat once on this PC before first start.
  set "PY=%PORTABLE_PY%"
) else (
  echo [ERROR] No Python found. Put portable_python\ here and run setup_env.bat
  pause
  exit /b 1
)

echo Using: %PY%
"%PY%" -c "import flask,Bio,pandas,primer3" >nul 2>&1
if errorlevel 1 (
  echo [ERROR] Dependencies missing. Run setup_env.bat first.
  pause
  exit /b 1
)

rem Ensure local package pcrlib is importable (portable Python / isolated ._pth).
set "PYTHONPATH=%~dp0;%PYTHONPATH%"

echo Starting web UI at http://127.0.0.1:5000/
"%PY%" "%~dp0web_app.py"
if errorlevel 1 (
  echo.
  echo App exited with an error.
  pause
)
