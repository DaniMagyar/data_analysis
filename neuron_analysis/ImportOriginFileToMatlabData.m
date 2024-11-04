% ***********************************
% MIT License
% 
% Copyright (c) 2022 Ruben Millan-Solsona
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
%%***********************************

function outData=ImportOriginFileToMatlabData(filepath)
%% Function that returns a structure with the books and sheets of the specified OriginPro document
%   If you specify a correct file path the function will return a
%   structure with the data from the Origin. Otherwise there would be
%   a dialog window for you to specify the *.opju or *opj file you want to load.
%%-----------------------------------------------------------------------------------

    %% Default input arguments
    if nargin<1 || ~ischar(filepath) || ~exist(filepath,'file')
        filepath = []; % Default is to open a file selection dialog
    end
    %% Display warning dialog
    choice = questdlg(sprintf('%s\n\n%s%s\n',...
    'Any currently open Origin Projects will be saved and closed in order to continue.',...
    'If an open Origin project does not have a filepath specified already the project will be saved ',...
    'to a new file in your Origin User Files directory, entitled ''UNTITLED.opj'''),...
    'Warning!',...
    'OK','Cancel','OK');
    if strcmpi(choice,'Cancel') % If user pressed Cancel, return
        outData = []; % Reset data structure        
        return
    end
    
    %% Obtain Origin COM Server object. This will connect to an existing
    % instance of Origin, or create a new one if none exist
    try originObj = actxserver('Origin.ApplicationSI');

    catch err
        if strcmpi(err.identifier,'MATLAB:COM:InvalidProgid') % If an Origin server object was not found on the computer
            msgbox('Sorry, OriginPro must be installed on your computer.','importOrigin');
        end       
        outData = []; % Reset data structure             
        return
    end
    
    % Select Origin file
    if isempty(filepath)
        [FileName,PathName,chose] = uigetfile('*.opju','Open Existing Origin Project','OriginProject.opju');

        if chose == 0 % If user pressed Cancel
            outData = []; % Reset data structure   
            return
        else
            filepath = fullfile(PathName,FileName); % Full file path
        end
    end
    % Open waitbar, opening Origin takes a while
    hWaitbar = waitbar(0,'Opening Origin. Please wait...','Name','ImportOrigin');

    % Save and close existing projects
    invoke(originObj, 'Execute', 'Save');
    invoke(originObj, 'IsModified', 'false');
    
    % Load the specified project
    invoke(originObj, 'Load', strcat(filepath));
    set(originObj,'Visible',1);
   
    workbooksHandle = invoke(originObj,'WorksheetPages');
    nbooks = get(workbooksHandle,'Count'); % Number of workbooks in project
    outData.Handle=workbooksHandle;
    outData.nbooks=nbooks;
    for b = 0:nbooks-1 % Loop all workbooks

        workbookHandle = get(workbooksHandle,'Item',int32(b)); % Handle to workbook b
        workbookName = get(workbookHandle,'Name'); % Short name of workbook b
        workbookLongName = get(workbookHandle,'LongName'); % Long name of workbook b

        % Identify worksheets
        worksheetsHandle = get(workbookHandle,'layers'); % Handle to worksheets of workbook b
        nsheets = get(worksheetsHandle,'Count'); % Number of sheets in workbook b
        % Output structure load
        outData.books(b+1).Handle=workbookHandle;
        outData.books(b+1).Name=workbookName;
        outData.books(b+1).LongName=workbookLongName;       
        outData.books(b+1).nsheets=nsheets;
        
        % Update waitbar
        waitbar((b+1/nbooks),hWaitbar);
        
        for s = 0:nsheets-1 % Loop all worksheets
            worksheetHandle = get(worksheetsHandle,'Item',int32(s)); % Handle to sheet s
            worksheetName = get(worksheetHandle,'Name'); % Short name of sheet s
            worksheetLongName = get(worksheetHandle,'LongName'); % Long name of sheet s
            worksheetData = invoke(originObj,'GetWorksheet',sprintf('[%s]%s',workbookName,worksheetName)); % Get worksheet data

            outData.books(b+1).Sheets(s+1).Handle=worksheetHandle;
            outData.books(b+1).Sheets(s+1).Name=worksheetName;
            outData.books(b+1).Sheets(s+1).LongName=worksheetLongName;
            outData.books(b+1).Sheets(s+1).Data=worksheetData;
            
            % Identify columns
            columnsHandle = get(worksheetHandle,'Columns'); % Handle to columns of worksheet s
            ncolumns = get(columnsHandle,'Count'); % Number of columns in sheet s
            outData.books(b+1).Sheets(s+1).ncolumns=ncolumns;
                        
           
            for c = 0:ncolumns-1 % Loop all columns
                columnHandle = get(columnsHandle,'Item',int32(c)); % Handle to column c
                columnName = get(columnHandle,'Name'); % Name of column c
                columnType = get(columnHandle,'Type'); % Column type: 0 (Y), 3 (X) or 5 (Z)
                columnLongName = get(columnHandle,'LongName'); % Long name of column c
                columnsUnits = get(columnHandle,'Units'); % Units specified in column c (not actually used here)
                columnComment = get(columnHandle,'Comments'); % Comments specified in column c (not actually used here)

               outData.books(b+1).Sheets(s+1).Columns(c+1).Handle=columnHandle;
               outData.books(b+1).Sheets(s+1).Columns(c+1).Name=columnName;
               outData.books(b+1).Sheets(s+1).Columns(c+1).columnType=columnType;
               outData.books(b+1).Sheets(s+1).Columns(c+1).LongName=columnLongName;
               outData.books(b+1).Sheets(s+1).Columns(c+1).Units=columnsUnits;
               outData.books(b+1).Sheets(s+1).Columns(c+1).Comment=columnComment;


            end
        end
    end  
    close(hWaitbar);
    % Close Origin win
    invoke(originObj,'Exit');
end