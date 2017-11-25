function xls_protect_sheets(xlsfile,sheets,action,password)
% 
% Protect or unprotect selected worksheets in Excel file, with or without
% password
% 
% 
%USAGE
%-----
% xls_protect_sheets(xlsfile,sheets,action)
% xls_protect_sheets(xlsfile,sheets,action,password)
% 
% 
%INPUT
%-----
% - XLSFILE : name of the Excel file
% - SHEETS  : cell array with the worksheet names, or matrix with positive
%   integers, or 'all', to tell which worksheets are going to be protected
% - ACTION  : 'protect' or 'unprotect'
% - PASSWORD: password to protect or unprotect the worksheet
% 
% 
%OUTPUT
%------
% - XLSFILE will be edited
% 
% 
% See also XLS_CHECK_IF_OPEN
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2012-Dec-30, 05:33 pm


% Input
%==========================================================================
if exist(xlsfile,'file')~=2
    fprintf('%s not found.\n',xlsfile)
    return
end
if ~strcmp(action,'protect') && ~strcmp(action,'unprotect')
    disp('Unknown option for ACTION.')
    return
end
if nargin<4
    password = '';
end


% Close Excel file
%-----------------
tmp = xls_check_if_open(xlsfile,'close');
if tmp~=0 && tmp~=10
    fprintf('%s could not be closed.\n',xlsfile)
    return
end


% The full path is required for the command "Workbooks.Open" to work
% properly
%-------------------------------------------------------------------
if isempty(strfind(xlsfile,filesep))
    xlsfile = fullfile(pwd,xlsfile);
end


% Read Excel file
%==========================================================================
%[type,sheet_names] = xlsfinfo(xlsfile);      % get information returned by XLSINFO on the workbook
Excel      = actxserver('Excel.Application'); % open Excel as a COM Automation server
set(Excel,'Visible',0);                       % make the application invisible
 % or ExcelApp.Visible = 0;
set(Excel,'DisplayAlerts',0);                 % make Excel not display alerts (e.g., sound and confirmation)
 % or Excel.Application.DisplayAlerts = false; % or 0
Workbooks  = Excel.Workbooks;                 % get a handle to Excel's Workbooks
Workbook   = Workbooks.Open(xlsfile);         % open an Excel Workbook and activate it
Sheets     = Excel.ActiveWorkBook.Sheets;     % get the sheets in the active Workbook
num_sheets = Sheets.Count;                    % number of worksheets


% Get sheets to protect
%==========================================================================
if ischar(sheets)
    if strcmp(sheets,'all')
        sheetidx = 1:num_sheets;
    else
        sheets   = {sheets};
        sheetidx = 1;
    end
elseif iscell(sheets)
    sheetidx = 1:length(sheets);
elseif isnumeric(sheets) && ~isempty(sheets) && ...
        isvector(sheets)==1 && all(floor(sheets)==ceil(sheets))
    sheetidx = sheets;
    sheetidx(sheetidx<1) = []; % minimum sheet index is 1
    sheetidx(sheetidx>num_sheets) = []; % maximum sheet index is num_sheets
else
    disp('Invalid input for SHEETS.')
    xls_save_and_close
    return
end


% Protect sheets
%==========================================================================
try
    for ss=sheetidx
        if iscell(sheets)
            sheet = get(Sheets,'Item',sheets{ss});
        else
            sheet = get(Sheets,'Item',ss);
        end
        sheet.Activate;
        invoke(Excel.ActiveSheet,action,password);
    end
catch ME
    disp(ME.message)
    if ss==sheetidx(end)
        if iscell(sheets)
            tmp = sheets{ss};
        else
            tmp = Sheets.Item(ss).Name;
        end
        fprintf('"%s" could not be %sed.\n',tmp,action)
    else
        if iscell(sheets)
            tmp  = sheets{ss};
            tmp2 = sheets{sheetidx(end)};
        else
            tmp  = Sheets.Item(ss).Name;
            tmp2 = Sheets.Item(sheetidx(end)).Name;
        end
        fprintf('"%s" to "%s" could not be %sed.\n',tmp,tmp2,action)
    end
end

xls_save_and_close


% Save and close
%==========================================================================
    function xls_save_and_close
    Workbook.Save;   % save the workbook
    Workbooks.Close; % close the workbook
    Excel.Quit;      % quit Excel
    % or invoke(Excel,'Quit');
    delete(Excel);   % delete the handle to the ActiveX Object

    clear Excel Workbooks Workbook Sheets sheet
    
    end

end