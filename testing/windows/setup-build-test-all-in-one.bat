

set PREVCWD=%CD%
set PREVPATH=%PATH%

@if "%ARTIFACTS_DIR%"=="" (
	set ARTIFACTS_DIR=artifacts
	@echo ARTIFACTS_DIR was not defined, set to "artifacts" )

:: for TOOLING_DIR
call %~dp0\default-tooling-path

:: Clear the tooling if there is, set it up from scratch
rd /s /q %TOOLING_DIR%
call %~dp0\download-setup-requirements || goto :error

md %ARTIFACTS_DIR%

:: Build for 32bit, test, and get wheel
set WHEEL_TO_TEST=%REPO_DIR%\testing\sandbox\python_package\dist\eden_simulator-%VERSION_PY%-py3-none-win32.whl
call %~dp0\build-eden-i686.bat || goto :error
call %~dp0\run-tests-on-wheel.bat %WHEEL_TO_TEST% || goto :error

copy %WHEEL_TO_TEST% %ARTIFACTS_DIR% || goto :error

:: Prepare for a new build, reset PATH and cwd
set PATH=%PATH%
cd %PREVCWD%

:: Build for 64bit, test,and get wheel
set WHEEL_TO_TEST=%REPO_DIR%\testing\sandbox\python_package\dist\eden_simulator-%VERSION_PY%-py3-none-win_amd64.whl
call %~dp0\build-eden-amd64.bat || goto :error
call %~dp0\run-tests-on-wheel.bat %WHEEL_TO_TEST% || goto :error

copy %WHEEL_TO_TEST% %ARTIFACTS_DIR% || goto :error%ARTIFACTS_DIR% || goto :error


:: done !
exit /b 0

:error
@echo ERROR : Last command failed with error code %errorlevel%. Exiting.
@exit /b %errorlevel%
