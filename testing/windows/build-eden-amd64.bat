
call %~dp0\default-repo-path.bat || goto :error

:: add tooling to PATH
call %~dp0\setpath-amd64 || goto :error

:: determine VERSION to tag
call %~dp0\get-version || goto :error
set VERSION_PY=%VERSION%

@echo on
pushd %REPO_DIR%
make clean || goto :error
make -j%NUMBER_OF_PROCESSORS% eden wheel CFLAGS_extra="-static" BUILD=release BUILD_STAMP=%VERSION% WHEEL_PLAT=win_amd64 WHEEL_VERSION=%VERSION_PY% || goto :error

popd

@goto :EOF

:error
@echo ERROR : Last command failed with error code %errorlevel%. Exiting.
@exit /b %errorlevel%
