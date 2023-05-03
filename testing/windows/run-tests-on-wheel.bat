
@if "%1"=="" (
	@echo Command needs a command-line argument for wheel to test. Exiting.
	@exit /b 2 )

set TEST_VENV_PATH=venv-eden-test

:: Use a clean venv
rmdir /q /s "%TEST_VENV_PATH%"
python3 -m virtualenv "%TEST_VENV_PATH%" || goto :error
call "%TEST_VENV_PATH%\Scripts\activate.bat" || goto :error_venv

python3 -m pip install wheel || goto :error_venv
python3 -m pip install %1 || goto :error_venv

call %~dp0\run-basic-tests.bat || goto :error_venv

call "%TEST_VENV_PATH%\Scripts\deactivate.bat" || goto :error

exit /b 0

:error_venv
call "%TEST_VENV_PATH%\Scripts\deactivate.bat"
:error
@echo ERROR : Last command failed with error code %errorlevel%. Exiting.
@exit /b %errorlevel%
