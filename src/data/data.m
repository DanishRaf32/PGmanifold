function varargout = data(D, action,varargin) %#ok<INUSL>
%DATA Create and manage the dataset table view
%
%  CB = DATA(OBJ, 'get_callbacks') returns a structure of function handles
%  to the functions that can be used to create and manage the dataset table
%  view.

%  Copyright 2000-2016 The MathWorks, Inc. and Ford Global Technologies, Inc.



if nargin > 1
    action = convertStringsToChars(action);
end

if nargin>1 && ischar(action)
    switch lower(action)
        case 'get_callbacks'
            varargout{1} = i_GetCallbacks;
        case 'createpage'
            [varargout{1:2}] = createPage(varargin{:});   
    end
end

%-----------------------------------------------------------------------
function cb = i_GetCallbacks
%-----------------------------------------------------------------------
cb.View = @i_View;
cb.Copy = @i_Copy;
cb.Show = @i_Show;
cb.Paste = @i_Paste;
cb.ExprClick = @i_ExprClick;
cb.Draw = @i_DrawDataPage;
cb.Selection = @i_Selection;
cb.ShowIndex = @i_ShowIndex;


function [d,newview]=createPage(d)

cb = i_GetCallbacks;

% -------- data page
newview = struct('ID','data',...
    'label','&Data',...
    'icon','cgdsdatabut.bmp',...
    'tooltip','View Data',...
    'drawcb',cb.Draw,...
    'layout',1,...
    'view',cb.View,...
    'show',cb.Show,...
    'copy',cb.Copy,...
    'paste',cb.Paste,...
    'tpmenu',d.Handles.fm.menu,...
    'tpselection',cb.Selection,...
    'showindex',cb.ShowIndex,...
    'tplist',0,...
    'bmlist',1,...
    'bmmenu',d.Handles.fm.menu,...
    'tlenable',d.Handles.tm.FactorBits,...
    'bmclick',cb.ExprClick);

%------------------------------------------------------------------
function [d,lyt] = i_DrawDataPage(d)
%------------------------------------------------------------------
Handles = d.Handles;
fig = Handles.Figure;

% Context menu for row header clicks
m = uicontextmenu('Parent', fig);
uimenu(m , 'Label' , '&Insert above',...
    'Callback' , {@cb_InsertRow, 0});
uimenu(m , 'Label' , 'Insert &below',...
    'Callback' , {@cb_InsertRow, 1});
uimenu(m , 'Label' , '&Delete',...
    'Callback' , @cb_RemoveRow);
Handles.TableRowMenu = m;

% Table
Handles.Table = cgdatasetgui.Table('parent', Handles.ViewParent,...
    'editable', true, ...
    'visible', 'off', ...
    'RowHeaderContextMenu', m);
Handles.TableEditListener = mbcgui.util.listener(Handles.Table, 'DataChanged', @cb_EditClick);
Handles.TableSelListener = mbcgui.util.listener(Handles.Table, 'SelectionChanged', @cb_SelChange);
Handles.TableColumnRCListener = mbcgui.util.listener(Handles.Table, 'ColumnPopup', @cb_OpenColContext);

lyt = Handles.Table;
d.Handles = Handles;


%------------------------------------------------------------------
function d = i_Show(d)
%------------------------------------------------------------------
page = d.ViewInfo(d.currentviewinfo);

d.Handles.FactorList.SelectionCallback = page.tpclick;
d.Handles.ExprList.SelectionCallback = page.bmclick;

d.Handles.FactorList.ContextMenu = page.tpmenu;
d.Handles.ExprList.ContextMenu = page.bmmenu;

configureExprList(d.Handles.ExprList);
d.Handles.fm.BmVis = d.Handles.fm.FactorVisBm;


%------------------------------------------------------------------
function d = i_View(d,sel_name)
%------------------------------------------------------------------
% Completely redraw the table
d.Handles.Table.setDataset(d.pD);
set(d.Handles.TopCard,'currentcard',d.currentcard);

if nargin<2
    sel_name = -1;
end
d = pr_RefreshExprList(d,sel_name,'cage');

d = i_ExprClick(d, d.Handles.FactorList);


%------------------------------------------------------------------
function cb_InsertRow(~, ~, pos)
%------------------------------------------------------------------
d = pr_GetViewData;
t = d.Handles.Table;

if nargin<3
    pos = 1;
end

row = double(t.getSelectedRows);
nRows = length(row);
if nRows
    if pos==1
        row = max(row);
    else
        row = min(row)-1;
    end
else
    nRows = 2;
    row = 1;
end

% Add new row to underlying dataset
d.pD.info = d.pD.AddRow([], row, nRows);

% Perform quick refresh of table cell data
t.refresh;


%------------------------------------------------------------------
function cb_RemoveRow(~, ~)
%------------------------------------------------------------------
d = pr_GetViewData;
t = d.Handles.Table;

row = double(t.getSelectedRows);
if ~isempty(row)
    OpPt = d.pD.info;
    
    % Remove data from underlying dataset
    data = get(OpPt, 'Data');
    data(row,:) = [];
    OpPt = set(OpPt, 'Data',data);
    
    % Convert grids into a block of data
    OpPt = convertGridToBlock(OpPt);
    d.pD.info = OpPt;
    
    % Perform quick refresh of table cell data
    t.refresh;
end

%------------------------------------------------------------------
function cb_EditClick(~,~)
%------------------------------------------------------------------
% Update the toolbar when the dataset is edited
d = pr_GetViewData;
pr_EnableToolbar(d);

%------------------------------------------------------------------
function cb_SelChange(~, ~)
%------------------------------------------------------------------
% Deselect the items in the project expression list
d = pr_GetViewData;
deselect(d.Handles.ExprList);


%------------------------------------------------------------------
function cb_OpenColContext(~, ~)
%------------------------------------------------------------------
% This is general code for displaying context menus
d = pr_GetViewData;
pl=get(0,'PointerLocation');
fp=get(d.Handles.Figure,'Position');
fp=fp(1:2);

page = d.ViewInfo(d.currentviewinfo);
list = 'top';

Menu = page.tpmenu;
MenuCB = get(Menu,'Callback');

top_sel = feval(page.tpselection,d);
feval(MenuCB,d,list,top_sel);

set(Menu,'Visible','off');
pos = pl-fp;
set(Menu,'Position',pos);
set(Menu,'Visible','on');

%-----------------------------------------------------------------------
function sel = i_Selection(d)
%-----------------------------------------------------------------------
% Only report a selection when we are trying to open the column context
% menu
sel = d.Handles.Table.getSelectedDatasetColumns;
sel = d.Exprs.shown_factors(sel);


%-----------------------------------------------------------------------
function d = i_ShowIndex(d,index)
%-----------------------------------------------------------------------
% Highlight the thing given in index (prob. an error column).
d.Handles.Table.selectDatasetColumns(index);

%-----------------------------------------------------------------------
function d = i_ExprClick(d,~)
%-----------------------------------------------------------------------
% Click on exprlist
hList = d.Handles.ExprList;
index = hList.Selected;
if ~isempty(index)
    % Select appropriate column of table.
    index = index(index<=length(d.Exprs.factor_index));
    % Convert to true dataset index ordering
    oppt_index  = d.Exprs.factor_index(index);
    oppt_index = oppt_index(oppt_index>0);
    if isempty(oppt_index)
        d.Handles.Table.clearSelection;
    else
        d.Handles.Table.selectDatasetColumns(oppt_index);
    end
else
    d.Handles.Table.clearSelection;
end


%------------------------------------------------------------------
function data = i_Copy(d)
%------------------------------------------------------------------
% Copying can be delegated entirely to the table object
if ~isempty(getSelectedColumns(d.Handles.Table))
    d.Handles.Table.copy;
end
data = [];

%------------------------------------------------------------------
function i_Paste(d,~)
%------------------------------------------------------------------
% Pasting can be delegated entirely to the table object
d.Handles.Table.paste;
