@echo off
setlocal EnableExtensions
cd /d "%~dp0"

set "PY_PORTABLE=%~dp0portable_python\python.exe"
set "VENV_DIR=%~dp0.venv"
set "PTH=%~dp0portable_python\python314._pth"
set "GETPIP=%~dp0portable_python\get-pip.py"

echo ============================================
echo  Setup portable environment for In silico PCR
echo ============================================

if not exist "%PY_PORTABLE%" (
  echo [ERROR] portable_python\python.exe not found.
  echo Put Windows embeddable Python into portable_python\
  exit /b 1
)

echo [1/4] Enable site-packages in portable Python...
powershell -NoProfile -Command ^
  "$p='%~dp0portable_python\python314._pth'; if (Test-Path $p) { $c=Get-Content $p -Raw; if ($c -notmatch '(?m)^import site\s*$') { $c=$c -replace '(?m)^#\s*import site\s*$','import site'; if ($c -notmatch '(?m)^import site\s*$') { $c=$c.TrimEnd()+\"`r`nimport site`r`n\" }; Set-Content -Path $p -Value $c -NoNewline } }"

echo [2/4] Ensure pip...
"%PY_PORTABLE%" -m pip --version >nul 2>&1
if errorlevel 1 (
  echo pip not found, downloading get-pip.py ...
  powershell -NoProfile -Command ^
    "Invoke-WebRequest -Uri 'https://bootstrap.pypa.io/get-pip.py' -OutFile '%GETPIP%'"
  if not exist "%GETPIP%" (
    echo [ERROR] Failed to download get-pip.py. Check internet access.
    exit /b 1
  )
  "%PY_PORTABLE%" "%GETPIP%"
  if errorlevel 1 (
    echo [ERROR] get-pip failed.
    exit /b 1
  )
)

echo [3/4] Create virtual environment .venv ...
if exist "%VENV_DIR%\Scripts\python.exe" (
  echo .venv already exists, reuse it.
) else (
  "%PY_PORTABLE%" -m venv "%VENV_DIR%"
  if errorlevel 1 (
    echo [WARN] python -m venv failed. Falling back to installing into portable_python.
    set "VENV_PY=%PY_PORTABLE%"
    goto INSTALL
  )
)

set "VENV_PY=%VENV_DIR%\Scripts\python.exe"

:INSTALL
echo [4/4] Install requirements into environment...
"%VENV_PY%" -m pip install --upgrade pip
"%VENV_PY%" -m pip install -r "%~dp0requirements.txt"
if errorlevel 1 (
  echo.
  echo [ERROR] Dependency install failed.
  echo Portable Python is 3.14 — some wheels ^(primer3-py / biopython^) may be missing.
  echo Prefer portable Python 3.11 or 3.12 if this happens.
  exit /b 1
)

echo.
echo ============================================
echo  Setup OK
echo  Run: run_web_app.bat
echo ============================================
exit /b 0
