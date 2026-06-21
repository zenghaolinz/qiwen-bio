@echo off
setlocal

cd /d "%~dp0"
title Qiwen Bio

set "VENV_PYTHON=%CD%\.venv\Scripts\python.exe"

if not exist "%VENV_PYTHON%" (
    echo [Qiwen Bio] Creating Python virtual environment...
    where python >nul 2>nul
    if errorlevel 1 (
        echo [Qiwen Bio] Python was not found in PATH.
        echo Install Python 3.11 or newer, then run this script again.
        goto :error
    )
    python -m venv .venv
    if errorlevel 1 goto :error
)

"%VENV_PYTHON%" -c "import fastapi, uvicorn, qiwen_bio" >nul 2>nul
if errorlevel 1 (
    echo [Qiwen Bio] Installing project dependencies...
    "%VENV_PYTHON%" -m pip install -e .
    if errorlevel 1 goto :error
)

echo [Qiwen Bio] Starting http://127.0.0.1:8000
echo [Qiwen Bio] Press Ctrl+C to stop the server.
"%VENV_PYTHON%" -m uvicorn qiwen_bio.api:app --host 127.0.0.1 --port 8000 %*
if errorlevel 1 goto :error
goto :eof

:error
echo.
echo [Qiwen Bio] Startup failed.
pause
exit /b 1
