function aviobj = safe_slover_save(obj, opt_slvr_save)
% 
% SAFE_SLOVER_SAVE saves the slice view images.
% 
% 
%USAGE
%-----
% aviobj = safe_slover_save(obj, opt_slvr_save)
% 
% 
%INPUT
%-----
% - OBJ: slover object
% - OPT_SLVR_SAVE: structure with the following fields:
%   - SAVE_SLC: 1 (create image with slices) or 0 (don't create)
%     - If SAVE_SLC=1:
%       SLC_OUTPUT: output file name (full path, including extension, e.g.,
%       bmp, jpg, png, tif)
%       WRITE_TXT: 1 (write text above the first slice) or 0 (don't write)
%       - If WRITE_TXT=1:
%         TXT: text that should go above the first slice
%   - SAVE_MOVIE: 1 (create AVI file) or 0 (don't create)
%     By editing the script it is possible to create animated GIFs
%     - If SAVE_MOVIE=1:
%       DELAY     : time between the task onset and beginning of HRF
%       CONDNAME  : condition name for that frame (string [can be empty])
%       MOVIE_FILE: output movie file name (full path)
%       FRAMEDUR  : frames duration [s]
% 
% 
%OUTPUT
%------
% - If SAVE_SLC=1:
%   - SLC_OUTPUT
% - If SAVE_MOVIE=1:
%   - MOVIE_FILE
%   - AVIOBJ: AVI object created with function "avifile"
% 
% 
% Based on slover.m (slice overlay), pr_basic_ui.m
% 
% See also SAFE_SLOVER_CREATE
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2014-Jul-14


% Read options
%==========================================================================
save_slc   = opt_slvr_save.save_slc;
save_movie = opt_slvr_save.save_movie;
if save_slc
    output    = opt_slvr_save.slc_output;
    write_txt = opt_slvr_save.write_txt;
    if write_txt
        txt = opt_slvr_save.txt;
    end
end
if save_movie
    delay      = opt_slvr_save.delay;
    condname   = char(opt_slvr_save.cond); % must be a string
    movie_file = opt_slvr_save.movie_file;
    framedur   = opt_slvr_save.framedur;
end
visible = 'off'; % figure visibility: 'on' or 'off'


% Use SPM figure window
%--------------------------------------------------------------------------
obj.figure = spm_figure('GetWin', 'Graphics');
bg_color = [0 0 0];
set(obj.figure, 'Color'  , bg_color) % change background color to black
set(obj.figure, 'Visible', visible ) % figure visibility


% Save figure
%==========================================================================
if save_slc
    
    
    % Display slices
    %----------------------------------------------------------------------
    %if write_txt
    %    obj.area.position = [0 0 1 0.9600]; % adjust position to fit text above
        % default: [0 0 1 0.9200]
    %end
    obj = paint(obj);
    set(obj.figure, 'Visible', visible) % figure visibility
    
    
    % Change font size in the colorbar
    %----------------------------------------------------------------------
    %set(findobj('Parent',obj.figure,'Tag','cbar'),...
    %    'FontSize',.1,'FontWeight','bold')
    
    
    % Guarantee that an old file is not overwritten
    %----------------------------------------------------------------------
    tmp = 1;
    if exist(output,'file')==2
        ext    = output(end-2:end);
        output = sprintf('%s_%.2d.%s',output(1:end-4),tmp,ext);
        while exist(output,'file')==2
            tmp = tmp + 1;
            output(end-5:end-4) = sprintf('%.2d',tmp);
        end
    % else % do nothing
    end
    
    
    % Write text above the first slice
    %----------------------------------------------------------------------
    if write_txt
        write_txt_paint(obj.figure, txt, 'left')
    end
    
    
    % Save figure
    %----------------------------------------------------------------------
    %set(obj.figure,'PaperPositionMode','auto') % this line corrects the fact that sometimes the output figure is cut
    drawnow; pause(.1)
    %saveas(obj.figure,output,output(end-2:end))
    
    % Size of the SPM Graphics window
    %Rect = [515 15 600 865]; Rect = Rect.*spm('WinScale'); Rect = Rect(3:4);
    fig_print(obj.figure, output, [1080 1557], 1, 1, bg_color)
    
	% Print the output file name
	%----------------------------------------------------------------------
	% Useful too see where the files are being saved
    %fprintf('%s\n',output)
    
end


% Create movie
%==========================================================================
if save_movie
    
    
    % Do the display
    %----------------------------------------------------------------------
    %obj.area.position = [0 0 1 0.9600]; % adjust position to fit text above
    % default: [0 0 1 0.9200]
    obj = paint(obj);
    set(obj.figure,'Visible',visible) % figure visibility
    
    
    % Change font size in the colorbar
    %----------------------------------------------------------------------
    %set(findobj('Parent',obj.figure,'Tag','cbar'),...
    %    'FontSize',.1,'FontWeight','bold')
    
    
    % Write time above the 1st slice and condition above the 3rd slice (if
    % possible)
    %----------------------------------------------------------------------
    write_txt_paint(obj.figure,{sprintf('%5.1f s',delay),condname},'center')
    
    
    % Convert RGB to indexed image and get the colormap
    %----------------------------------------------------------------------
    movegui(obj.figure) % in case figure goes to the 2nd monitor
    drawnow; pause(.1)
    fig = getframe(obj.figure);
    % For +ve (-ve) BOLD, blue (red) gets gray. The command below works
    %[img,map] = rgb2ind(fig.cdata,96,'nodither'); % for animated GIF
    % 96 is the minimum number of colors needed
    % Images are practically the same with dither, but the files are larger
    
    
    % Create or update movie file
    %----------------------------------------------------------------------
    
    % Animated GIF
%     try
%         % returns error if it is the first frame
%         imwrite(img,map,movie_file,'DelayTime',framedur,'WriteMode','append')
%     catch
%         imwrite(img,map,movie_file,'DelayTime',framedur,'LoopCount',inf) % first AVI frame
%     end
    
    % AVI
    if exist(movie_file,'file')==2 % not the first frame
        %aviobj = addframe(opt_slvr_save.aviobj,obj.figure);
        aviobj = addframe(opt_slvr_save.aviobj,fig.cdata);
    else % first frame
        aviobj = avifile(movie_file,...
            'compression','None',...
             'fps',1/framedur);
        % 'colormap': m-by-3 matrix, m<=256 (236 for Indeo compression)
        % 'compression' (Windows): 'MSVC'; 'RLE'; 'Cinepak' (32-bit);
        %   'Indeo3' or 'Indeo5' (32-bit Windows XP); 'None'
        % 'keyframe'
        % 'quality': 0-100 (only when there is compression)
        %aviobj = addframe(aviobj,obj.figure);
        aviobj = addframe(aviobj,fig.cdata);
    end
    
end


% Change background color back to white
%--------------------------------------------------------------------------
set(obj.figure, 'Color', [1,1,1])

end


%==========================================================================


function write_txt_paint(fig,txt,pos)
% Write the string TXT at the top of the image at position POS ('left',
% 'center', 'right'). When TXT is a 2x1 cell array, cell 1 will be written
% over the 1st slice, and cell 2, over the 3rd slice if possible (or 2nd).

slc = findobj(fig,...
    'Type', 'axes',...
    'Tag' , 'slice overlay panel');
% axes(slc(end)) % select the first slice

% Get figure height (pixels)
%---------------------------
oldunits = get(fig, 'Units');
set(fig,'Units', 'pixels')
height = get(fig, 'Position');
height = height(4);

% Set font size
%--------------
oldslcunits = get(slc(end), 'Units');
set(slc(end), 'Units', 'pixels')
font = get(slc(end), 'Position');
font = font(2) + font(4); % top position of the slice (pixels)
font = height - font;     % size of the top margin (pixels)
font = font/1.5; % this seems to return a good font size

% Adjust font size
%-----------------
height = get(slc(end), 'Position');
height = height(4); % slice height (pixels)
if font/height>0.25 % correct font size if it is too big
    font = 0.25*height;
end

% Check if cells can be written over different slices
%----------------------------------------------------
if iscellstr(txt)
    slc_pos = [0 0];
    tmp     = get(slc(end), 'Position');
    if length(slc)>2
        tmp2 = get(slc(end-2), 'Position');
        if abs(tmp(2)-tmp2(2))<10^(-8) % same vertical position
            slc_pos(2) = 2;
        else
            tmp2 = get(slc(end-1), 'Position');
            if abs(tmp(2)-tmp2(2))<10^(-8) % same vertical position
                slc_pos(2) = 1;
            else
                %slc_pos(2) = 0;
                txt{1} = [txt{1} '  ' txt{2}];
                txt(2) = [];
            end
        end
    else
        %slc_pos(2) = 0;
        txt{1} = [txt{1} '  ' txt{2}];
        txt(2) = [];
    end
    set(slc(end-slc_pos(2)), 'Units', 'pixels')
else
    txt     = {txt};
    slc_pos = 0;
end

% Write text
%-----------
for tt=1:length(txt)
    
    % Title
    title(slc(end-slc_pos(tt)), txt{tt},...
        'Color'     ,'w',...
        'FontUnits' ,'pixels',...
        'FontSize'  ,font,...
        'FontWeight','bold',...
        'HorizontalAlignment',pos)
    
    % Reset properties
    set(slc(end-slc_pos(tt)),'Units',oldslcunits)
    
end
    
% Reset properties
%-----------------
set(fig,'Units',oldunits)

end