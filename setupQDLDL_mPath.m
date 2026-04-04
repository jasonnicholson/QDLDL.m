function setupPath(shouldSetup)
% setupPath adds or removes the src folder from the MATLAB path
  arguments (Input)
    shouldSetup (1,1) logical = true;
  end

  repoRoot = fileparts(mfilename("fullpath"));
  src = fullfile(repoRoot, "src/");

  if shouldSetup
    addpath(src);
  else
    rmpath(src);
  end
end