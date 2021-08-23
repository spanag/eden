
:: for TOOLING_DIR
call %~dp0\default-tooling-path

:: Clear the tooling if there is, set it up from scratch
rd /s /q %TOOLING_DIR%
call %~dp0\download-setup-requirements || goto :error

call %~dp0\build-test-all-in-one

:: done !
exit /b 0

:error
@echo ERROR : Last command failed with error code %errorlevel%. Exiting.
@exit /b %errorlevel%
