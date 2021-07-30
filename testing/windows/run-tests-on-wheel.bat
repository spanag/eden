
@if "%1"=="" (
	@echo Command needs a command-line argument for wheel to test. Exiting.
	@exit /b 2 )

python3 -m pip uninstall -y eden-simulator || goto :error
python3 -m pip install %1 || goto :error
call %~dp0\run-basic-tests.bat || goto :error

exit /b 0

:error
@echo ERROR : Last command failed with error code %errorlevel%. Exiting.
@exit /b %errorlevel%
