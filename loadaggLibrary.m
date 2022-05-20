%
% Load fortran library. If it is already loaded, it returns the library
% name.
%
function sLibname = loadaggLibrary

path = fileparts(mfilename('fullpath'));

    sLibname = 'aggcal';
    sExtension = '.dll';

    if ~libisloaded(sLibname)
        [notfound,warnings] = loadlibrary(...
            strcat(path,'\',sLibname,sExtension), ...
            strcat(path,'\aggcal_wrap.h'));
    end
end