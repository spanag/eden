
call %~dp0\default-repo-path.bat || goto :error

:: replace \ with / as seen on https://stackoverflow.com/a/23544993
set "REPO_DIR_SLASH=%REPO_DIR:\=/%"
:: echo "%REPO_DIR_SLASH%"

python3 -c "import eden_simulator; eden_simulator.runEden('%REPO_DIR_SLASH%/examples/LEMS_NML2_Ex25_MultiComp.xml')" || goto :error

exit /b 0

:error
@echo ERROR : Last command failed with error code %errorlevel%. Exiting.
@exit /b %errorlevel%
