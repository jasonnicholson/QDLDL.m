function results = run_all_tests
%RUN_ALL_TESTS Execute all MATLAB unit tests for the MATLAB port.
    thisDir = fileparts(mfilename('fullpath'));
    addpath(fileparts(thisDir));
    addpath(thisDir);
    results = runtests(thisDir);
end
