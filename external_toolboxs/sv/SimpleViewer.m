function SimpleViewer(img, CLim, roi, txt, tileSize, XLim, YLim, cmap, prefs)
% SimpleViewer  A simple keyboard controlled movie viewer
% 
% ----- Syntax -----------------------------------------------------------------
%  SimpleViewer(img [, CLim, roi, txt, tileSize, prefs])
%  SimpleViewer(fileName)
%
%  img      - can be either:
%             - a 5d matrix: img(:,:,slice,time,subplot)
%             - a 1d cell array: img{time}
%             - a 3d nested cell array: img{subplot}{slice}{time}
%  CLim     - colormap limits (analogous to CLim in imagesc)
%  txt      - image titles (format: txt{subplot}{slice}{time})
%             It is permissible to have only one title for all times within a
%             given slice and subplot.
%  roi      - roi vertices (format: roi{subplot}{slice}{time}{roi})
%             Each ROI is defined by vertices in a N-by-2 arrays
%             [X1 Y1;...;XN YN].
%  tileSize - tiling size for subplots (format: [nRows nCols])
%  XLim     - x-axis limits, either as a 1x2 array or a cell array of 1x2 arrays
%             specifying XLims for each subplot
%  YLim     - y-axis limits with the same format as XLim
%  cmap     - colormap, either a single Nx3 array or a cell array of Nx3 arrays
%             specifying colormaps for each subplot
%  prefs    - struct with entries that override default preferences
%
%  All arguments aside from img are optional and may be set as an empty matrix
%  in order to set arguments further down the list.  prefs may also be passed as
%  the last argument, even if it is not in the 6th position.  SimpleViewer may
%  be called with no arguments to display the calling syntax and a list of
%  preferences.
%
%  SimpleViewer may be called by passing the file name of a DICOM, NIfTI, any
%  image format supported by imread (e.g. bmp, jpg, gif, ...), or any video
%  format supported by VideoReader (e.g. avi, mpg, ...).  Wildcards can be used
%  to load multiple images/NIfTI files.
%
%  SimpleViewer may also be passed a SimpleViewer workspace (.svmat or .mat)
%  file created by using 'shift + s' within a previous SimpleViewer instance.
%
% ----- Overview ---------------------------------------------------------------
%  The most basic functionality of SimpleViewer is provide easy visualization of
%  an image series, often one or more sets of time series.  If these time series
%  are of consistent enough dimensions that they can be represented as an <= 5
%  dimensional matrix, SimpleViewer can simply be called by passing this matrix
%  in the format of (:,:,slice,time,subplot), where the first two dimensions are
%  y and x respectively (i.e. standard row/column).  If a 3D matrix is passed,
%  the third dimension is assumed to be time, unless the autoSlicesToCine pref
%  option is set to false.  The left/right arrow keys and mouse scroll may be
%  used to advance/reverse in time, and the up/down arrow keys are used to
%  change slice.  The time/slice changes are applied to all subplots.  The
%  spacebar may be used to return to a ("home") frame (default: 1).  'ctrl +
%  space' may be used to start/stop playing the data in a loop in the time
%  dimension.
%
%  The displayed contrast of the images may adjusted using the window/level tool
%  ('w'), and the zoom of the images can be adjusted using the pan/zoom tool
%  ('z'), as detailed below.  The colormap can be cycled using 'm' and a
%  colorbar displayed using 'b'.  The overall window size can be adjusted by
%  resizing the window or by using the +/- keys.
%
%  The second major feature of SimpleViewer is the tools to draw and manipulate
%  regions of interest (ROIs).  Arbitrary polygons can be drawn and then
%  manipulated by moving their center or vertices or scaling them
%  larger/smaller.  A nudge tool can also be used to move vertices with a
%  circular brush.  ROI editing can be done on subplots independently or linked
%  between subplots or copy/pasted between subplots or SimpleViewer windows.
%
%  The third aspect of SimpleViewer is to provide some basic analysis of image
%  data.  Information about an ROI such as its circumference, area, as well as
%  the mean and standard deviation of pixels contained within it can be toggled
%  using the 'i' key.  The 'g' key can be used to open up a secondary window
%  with further information and contains tables of the mean and standard
%  deviation and a plot of the mean +/- standard deviation signal intensity as a
%  function of time.  If a line ROI is drawn, the profile across that line is
%  shown across time (similar to an echocardiography 'm-mode').
%
%  Finally, SimpleViewer data may be exported as images (e.g. .jpg, .png, .gif,
%  etc), animated GIFs, or movies (.avi or .mp4).  A SimpleViewer "workspace",
%  containing the image data, ROIs, window/level settings, etc. can be also be
%  saved in a .svmat format (which is a standard .mat file), and these
%  workspaces can reloaded by passing the filename to SimpleViewer.
%
% ----- Keyboard shortcuts -----------------------------------------------------
% General Controls:
%  'left'/'right'      - advance/reverse time
%  'scroll wheel'      - advance/reverse time
%  'up'/'down'         - advance/reverse slice
%  'keypad 0' or '/'   - go back to last time/slice
%  'space'             - return to home time frame (prefs.homeTime; default: 1)
%  'ctrl + space'      - toggle playing through time at 30 fps
%  'j'                 - toggle "jitter" mode, where a "yo-yo" movie is played 
%                        from 3 frames before the current time to 3 frames after
%  'esc'               - close window
%  '+' (or '=')        - increase zoom by 25%
%  '-' (or '_')        - decrease zoom by 25%
%  '\'                 - reset zoom to prefs.defaultScale (default: 3)
%                        - use 'shift' to reset zoom to 1
%  'm'                 - cycle colormaps in prefs.colormaps (default: gray/jet)
%  'alt' +             - select correspondingly numbered subplot
%  '1', '2', ... '9' 
%  ('shift') + 'tab'   - select next/previous subplot
%  'h'                 - hide figure text and SD bars on plot figure
%  'i'                 - cycle between ROI info display (none, ROI number, area,
%                        signal intensity, both)
%  'f'                 - cycle between fliplr, rot90(x,2), and flipud
%  'b'                 - toggle colorbar display
%  'shift + scroll'    - zoom in/out
%  'g'                 - toggle display of a secondary figure showing a plot of
%                        ROI signal intensities over time (and other info)
%  'q'                 - Unselect all ROIs and reset any modes (e.g. draw,
%                        pan/zoom, window/level, etc.)
%
% Export/Saving Controls:
%  'e'                 - export image/movie
%  'shift + s'         - save SimpleViewer workspace
%  'alt + c'           - copy image to clipboard
%
% Window/Level Mode Controls:
%  'w'                 - toggle window/level mode
%  'left click-drag'   - manually adjust window level (up/down for level,
%                        left/right for window). This window level is then
%                        applied to all subplots, unless the 'alt' is pressed.
%  'right click'       - automatic window level using the current subplot. This
%                        window level is then applied to all subplots.
%  'shift right click' - dynamic automatic leveling, where the automatic window
%                        level (as above) is re-applied every time the
%                        slice/time is changed
%  'double left click' - reset window level to default mode, where each subplot 
%                        has its own window level.
%  Note: if 'alt' is held during window level adjustment, changes are applied
%  only to the current subplot.
% 
% Pan/Zoom Mode Controls:
%  'z'                 - toggle pan/zoom mode
%  'left click-drag'   - adjust pan
%  'right click-drag'  - adjust zoom (up to zoom in; down to zoom out).  This
%                        applies to all subplots unless 'alt' is pressed (also
%                        true for pan)
%  'double left click' - reset pan/zoom to display whole image. Each subplot is
%                        has its own pan/zoom axes.
%  Note: if 'alt' is held during pan/zoom adjustment,  changes are applied only
%  to the current subplot.
% 
% ROI Mode Controls:
%  'd'                 - toggle ROI drawing mode
%  'left click'        - start a new ROI or add a vertex to an in-progress ROI
%  'right click'       - finish an in-progress ROI or create a single point ROI
%  'c'                 - copy ROI (use with 'shift' for global clipboard)
%  'v'                 - paste ROI (use with 'shift' for global clipboard)
%  'r'                 - replicate selected ROI across all times. If no ROIs are
%                        selected, all ROIs are replicated
%  's'                 - enter/exit "static ROI" mode, where the same ROI is
%                        used across all times within a slice
%  'l'                 - toggle linking ROI editing across subplots
%  't'                 - toggle show/hide all ROIs
%  'y'                 - toggle applying ROIs as masks to the data
%  'delete' OR         - if an ROI is currently being drawn, delete last vertex;
%  'backspace'           if a finished ROI is selected; delete selected ROI
%  shift + 'delete'    - delete all ROIs
%  or 'backspace'
%  '<' (or ',')        - shrink selected ROI(s) by prefs.roiScaleFactor
%  '>' (or '.')        - expand selected ROI(s) by prefs.roiScaleFactor
%                        - use the shift key to shrink/expand 5x faster
%  shift + arrows      - move ROIs by prefs.roiNudgeFactor pixels; default 0.5
%  '1', '2', ... '9'   - select correspondingly numbered ROI
%  '`' (or '~')        - select all ROIs
%
%  'p'                 - interpolate roi to have prefs.interpNPts points
%  'n'                 - nudge mode
%  'u'                 - undo last nudge
%  'o'                 - circle mode
%  '[' and ']'         - shrink/expand circle
%                        - use the shift key to shrink/expand 5x faster
%
% ----- Detailed Usage ---------------------------------------------------------
% Nested Cell Arrays
%  Image data can be supplied as a <=5 dimensional array if the size of images
%  is consistent over time, number of times is consistent across slices, and
%  number of slices is consistent across subplots.  However, there may be
%  instances in which this criteria is not satisfied, e.g. displaying images of
%  different sizes or displaying a time series of images along with a parametric
%  map.  In these cases, images must be supplied as a nested cell array, in
%  img{subplot}{slice}{time} format.  This allows flexibility in both image
%  size, time, and slice dimensions at the expense of increased complexity in
%  data formatting.  This format is required for ROIs and text inputs, as a
%  simple 5-dimensional array is not flexible enough to hold this data.  A
%  non-nested three-dimensional cell array may be considered in the future but
%  is currently not supported.
%
% Global Preferences
%  Preferences can also be set using the global variable "svPrefs", which takes
%  precendence over the defaultPrefs defined in this file, but are overridden by
%  the prefs argument when calling SimpleViewer.
%
% ROI Drawing
%  ROIs can be created while in the "draw" mode (toggled with 'd').  To draw a
%  single point ROI, simply right-click. To start a multi-point ROI, left-click
%  to draw the first vertex, and continue left-clicking to add additional
%  vertices.  The delete/backspace key can be used to delete the last vertex.
%  The ROI is completed with a right-click, which draws the last vertex and
%  finishes the ROI.  The ROI may also be completed with a shift right-click,
%  which draws the last vertex, but finishes the ROI as line segment.  Line
%  segments produce a spatial profile along the segment in the plot window (if
%  enabled).  They are represented internally as a list of vertices terminated
%  with NaNs.
%
% ROI Selection and Adjustment
%  Click on any part (center, vertices, edges) of an ROI to select it. The
%  currently selected ROI is indicated with dashed edge lines instead of solid
%  lines. The ROI may be moved by left-click dragging the center, and individual
%  vertices may also be moved in the same way. Shift left-click on a vertex to
%  delete it and shift left-click an edge to add a vertex.
%
% ROI Clipboard
%  The 'c' and 'v' keys can be used to copy and paste one or more ROIs to and
%  from a clipboard.  The clipboard is shared within a figure, and can be used
%  between subplots and frames (time or slice).  If the 'shift' modifier is used
%  with the 'c' and 'v' keys, ROIs are copied to/from a global clipboard shared
%  with all figures.
%
% Static ROI mode
%  When Static ROI is enabled, ROIs drawn are copied across all time frames
%  within the current slice.  In this mode, the time frame may also be changed
%  while drawing an ROI.  With Static ROI mode off, ROIs in each time frame are
%  individually adjusted.  A warning is displayed with Static ROI mode off
%  unless prefs.showStaticRoiStatus is set to false.
%
% Nudge Tool
%  Enter/exit the nudge tool by pressing 'n'.  In nudge mode, left click-hold
%  creates a circular "force field" around the mouse cursor.  With the button
%  held down, move the mouse around and any ROI vertex that touches the force
%  field is pushed away from the mouse cursor.  Right click-hold creates an
%  "inverse force field" that contains ROI vertices within its borders.  If an
%  ROI is selected, only that ROI is affected.  Otherwise, all ROIs are
%  affected.  The radius of the force fields is determined by the distance to
%  the nearest vertex when mouse button is depressed -- start with the mouse
%  cursor farther away from vertices to create a field with a larger radius.
%  The last change to an ROI with the nudge tool can be undone with the 'u' key,
%  where the "last change" is defined as the last left-click (drag), whether any
%  vertices were moved or not.
%
% Circle Tool
%  Enter/exit the circle tool by pressing 'o'.  A blue and red circle outline is
%  shown and follows the cursor.  The  '[' and ']' keys can be used to shrink 
%  and grow the circle respectively.  Left-click to place the ROI.  The default
%  radius and number points in the circle is controlled in the preferences via
%  'circleRadius' and 'circlePoints' respectively.
%
% Automatic Circle Tool
%  Circles can be automatically detected using the Hough transform via the
%  imfindcircles function in the Image Processing Toolbox. To use this within
%  SimpleViewer, window/level the image to generate good contrast between the
%  circles and its surroundings, then press shift+'o'.  Detected circles are
%  shown with a red outline.  Click on each circle in the order it should be
%  added and press enter when done.  If none of the detected circles are
%  appropriate, press enter to finish (and not esc).
%
% Interpolate ROI
%  Interpolate the selected ROI to have pref.interpNPts vertices, where the last
%  prefs.interpIgnoreLast edges of the ROI are ignored in the interpolation.  If
%  no ROIs are selected, all ROIs on the current frame are interpolated.  If
%  prefs.autoInterpolateRois is set to true, ROIs are automatically interpolated
%  as soon as it is finished being drawn and after nudging.
%
% Interpolate ROI (new)
%  - 'display' - Use prefs.interpAlgorithm to interpolate the displayed line
%                between ROI vertices
%  - 'finish'  - Upon finishing a new ROI, use prefs.interpAlgorithm to
%                interpolate it to pref.interpNPts points.  The ROI is not
%                further interpolated and lines are linear between vertices
%  - 'always'  - ROIs are interpolated when finished (i.e. 'finish' mode), and
%                also after each nudge. Note: This is the same as setting
%                prefs.autoInterpolateRois to true
%  - 'none'    - Do not automatically interpolate ROI vertices or displayed
%                lines between vertices
%
% Interpolate Across Time
%  For data with high temporal resolution or through-plane spatial resolution, 
%  it may be helpful to only trace ROIs on certain frames (e.g. every other 
%  frame) and linearly interpolate the ROIs on the remaining frames.  To do
%  this, trace on the desired frames, and press 'a' to interpolate the ROIs
%  to all other frames between the first and last frame with drawn ROIs.  It is
%  important to note that each vertex is linearly interpolated independently, so
%  ROIs of the same index number must have the same number of vertices in the
%  traced frames.  Additionally, each vertex should correspond to roughly the
%  same structure of interest in each traced frame.
%
% ROI Toggles
%  The display of ROIs can be toggled using the 't' key in order to facilitate
%  fine countour adjustments.  The 'y' key can be used to temporarily apply the
%  ROIs as a mask to visualize the pixels enclosed by each ROI.  If there are
%  2-4 ROIs, these are assumed to be short axis myocardial contours, i.e.
%  epicardium, endocardium, RV insertion point, and blood pool.
%
% ROI Raw Data Tables
%  The means and standard deviations within each ROI are displayed in tabular
%  format in the correspondingly named tabs of the ROI figure.  Values can be
%  copied from this table.
%
% ROI "M-mode"
%  For line segment ROIs (i.e. ROIs that are finished with a shift right-click)
%  where the image data has a time dimension, an "M-mode" like image is
%  shown in the "M-mode" tab of the ROI figure.  This image shows the signal
%  intensity profile along the line segment in the vertical dimension (top is
%  the beginning of the line segment) and along time in the horizontal dimension
%  (time increasing to the right).
%
% Saving Workspaces
%  SimpleViewer workspaces can be saved using the 's' key.  This saves all
%  variables necessary to restore a SimpleViewer session to a .svmat/.mat file.
%  Saved workspaces can be loaded by passing the .svmat/.mat file name when
%  calling SimpleViewer.  An .svmat file is identical in format to a standard
%  MATLAB .mat file, but may be automatically opened by SimpleViewer if an
%  appropriately defined opensvmat.m function is defined. e.g.:
%    function opensvmat(filename)
%      figure;
%      SimpleViewer(filename);
%    end
%
%  Note that similar functions can also be created for SimpleViewer supported
%  graphics and movie formats to enable drag-and-drop functionality, e.g.
%  opengif, openavi, etc functions.
%
% Export Image/movie
%  A SimpleViewer figure can be exported as a GIF, JPEG, or PNG image using the
%  'e' key.  A movie (in the time dimension) can be exported in animated GIF
%  format, AVI format ("Motion JPEG Video"), or mp4 format ("MPEG-4 Video").
%  Note that exporting avi/mp4 files from MATLAB versions prior to R2014b may
%  have dropped frames when encoded at ~30+ fps due to poor performance of the
%  VideoWriter internal function
%
% Plot Figure
%  If the preference 'suppressRoiPlot' is set to false, a new window is created
%  when an ROI is drawn. This window shows a plot of the signal intensity within
%  the ROI over the time dimension. If the preference 'showStdPlot' is set to
%  true, error bars are plotted at +/- one standard deviation and can be hidden
%  with the 'h' command. Alternatively, if a line segment ROI is drawn (i.e. an
%  ROI completed with a shift-click), a spatial profile along the line is drawn.
%  In both cases, the ROI may be adjusted and the plot will be updated in
%  real-time. Note that SimpleViewer does not properly support temporal and
%  spatial profiles simultaneously.
%
%  The mean and standard deviation of pixels within each ROI are displayed each
%  time frame in tabular format in the secondary plot figure.  Values can be
%  selected and copied to other spreadsheet programs.
%
% Home Times
%  The 'home times' is set using prefs.homeTime (default: 1).  Pressing spacebar
%  changes the displayed image to the home time.  If prefs.homeTime is an array,
%  pressing spacebar once changes to the first home time, again for the second
%  home time, etc.
%
% Preferences Editor
%  If the 'propertiesGUI' function is found in the path (available through the
%  MATLAB File Exchange), the F12 key can be used to edit the prefs struct for
%  an active SimpleViewer workspace
%
% ----- For SimpleViewer programmers -------------------------------------------
% ROI Transcription
%  By default, ROIs drawn on each subplot are independent.  In some cases, such
%  as with co-registered images or data with corresponding parametric maps, it
%  may be useful to draw on multiple subplots simultaneously.  This can be
%  controlled via the prefs.transcribeRoi with the following settings:
%   - false:       No transcribing (default)
%   - true:        Transcribe to all plots
%   - cell array:  Each cell is an array containing subplot indices.  Subplots
%                  in each array are linked and ROIs are copied between them
%  Transcribing between subplots can be toggled with the 'l' key if
%  prefs.transcribeRoi is not a cell array.
%
% ROI Colours
%  The order of ROI colours is controlled via prefs.roiColors, which may either
%  be a character array of pre-defined MATLAB colours (e.g. 'rgb' for red,
%  green, and blue ROIs), or an (n x 3) array of RGB colours (i.e. the format
%  for colourmaps).  ROI colours are cyclic, so the colours repeat when the
%  number of ROIs exceeds the number of defined colours.
%
% Black and White ROIs
%  ROIs may shown as black lines with white vertices, which may be useful when
%  tracing on coloured images.  This is controlled via prefs.blackWhiteRoi with
%  the following settings:
%   - false:  Coloured ROIs according to prefs.roiColors (default)
%   - true:   All ROIs are black and white
%   - array:  Entries are subplot indices where black and white ROIs are used
% 
% Interactive Usage
%  The 'enter' key is intentionally left unbound.  To use SimpleViewer
%  interactively to get drawn ROIs, window level, etc, use the following code:
%
%    waitfor(hFig, 'CurrentCharacter', char(13));  % char(13) is "enter"
%    set(hFig, 'currentch', char(1));  % reset 'CurrentCharacter'
%    % Code to run after user input
%
%  The second line may be omitted, but waitfor() will immediately return if
%  called again later in the code on the same figure.
%
% Custom Keyboard Callbacks
%  SimpleViewer functionality may be extended through the use of custom
%  functions (callbacks) that are called when keys are pressed.  These are
%  defined via prefs.kbCallbacks, which is a struct in which each field name is
%  a key and each field value is a function handle.  This function must take
%  the SimpleViewer figure handle as its only argument.  Note that default
%  keyboard callbacks may be disabled by using an empty field value.
%
% Programmatic SimpleViewer Control 
%  Many internal SimpleViewer functions are exposed via prefs.fcnHandles, which
%  is a struct where each field name is the SimpleViewer function name and the
%  field value is the function handle.  This can be used to control SimpleViewer
%  from within scripts (e.g. change time/slice, enter/exit modes, etc).
%
% Contrast Adjustment Modifiers
%  Contrast (window/level) adjustments on one subplot are automatically applied
%  to all other subplots unless the alt key is held.  This may not be desirable
%  if there the displayed images have very different contrasts and may be
%  controlled with prefs.linkWindowLevel with:
%   - true:       Contrast is adjusted for all subplots simultaneously (default)
%   - false:      Each subplot is contrast adjusted independently
%   - cell array: Each cell is an array containing subplot indices.  Subplots
%                 in each array have contrast adjustments linked
%  Some data have meaningful values, such as the zero point for
%  positive/negative data.  It may be desirable to adjust the contrast window
%  without adjusting level, such that the zero point maintains the same colour.
%  This can be enabled by setting prefs.disableLevelAdjust to true.
%
% Simplified 'img' Formats
%  In order make SimpleViewer more accessible for quick usage, 'img' supports
%  a 5D matrix or a simple 1D cell array instead of the "nested cell" format.
%  Additionally, if prefs.autoSlicesToCine is true, img data which has multiple
%  slices but only a single time is convered to img data with a single slice and
%  multiple times.  Any other arguments (CLim, roi, txt, ...) are optional or
%  can be provided as empty matrices ([]) if necessary.
%
% 'Extra' Data
%  It may sometimes be useful to store additional data in a SimpleViewer figure,
%  such as metadata for the images or data used by custom callback functions.
%  The 'extras' field with setappdata(hFig, ...) is saved when saving workspaces
%  and is restored when loading workspaces.
%
% Notes:
%  - Subplots are ordered left-to-right in the top row, then left-to-right in
%    the next row down, etc., as with the subplot(m,n,p) syntax.  'hAxes' (used
%    internally), is a 2D matrix with rows and columns REVERSED of the subplot
%    layout, such that linear indexing into 'hAxes' works as expected.
%  - SimpleViewer recycles the current figure, much like image() and imagesc().
%  - SimpleViewer uses the OpenGL renderer for speed. When saving in eps format,
%    use "set(gcf, 'Renderer', 'painters')" so that ROI objects are editable.
%    Note also that 'h' will set this mode.
%
% Tips:
%  - A nested cell array m{i}{j} can be converted into a 2d cell array n{i,j}:
%      n = cat(1, m{:});
%    This assumes that all m{i} have the same number of elements 'j' and that
%    m{i} is a single-row cell array.
%  - A SimpleViewer window and its ROI data plot window are programmed to close
%    each other if one closes.  If a window is closed by some other means (e.g.
%    if "close all" is used), then the other window won't close properly.  In
%    that case, use select the figure and use the command:
%        delete(gcf)
%    Alternatively, delete all open figures using the command:
%        delete(get(0,'Children'))
%
% -------------------------------------------------------------------------
% Bug tracker:
% - UpdateRoiFigure() should actually not draw deleted ROI data instead of just
%   having a line plot with no data
% - The profile and m-mode tabs on the ROI plot figure cannot be exported in
%   <R2014b.
%  - No range checking for transcribeRois could lead to crashes
% -------------------------------------------------------------------------

% M-lint warning suppression:
%#ok<*INUSD>   Input argument '' might be unused.
%#ok<*USENS>   Variable '' is used, but might be unset.
%#ok<*DSPS>    disp(sprintf(...)) can usually be replaced by fprintf(...)

% Revisions:
%  4.2.4   - New: Added colormaps as an argument, enabling support for separate
%            colormaps for each subplot and saving of colormaps with svmat files
%  4.2.3   - New: Add support for HDF5 files
%          - New: Use prefs.useDicomColormaps to update the colormap for each
%            subplot, if DICOM headers are stored in the extras field
%          - New: Slice/time navigation with the keyboard arrow keys apply only
%            to the current subplot if 'alt' is held
%  4.2.2   - New: If DICOM headers are stored in the extras field (e.g. through
%            ReadDicomDir), then the intersection between images displayed in
%            different subplots can be displayed with 'x'
%            New: If DICOM headers are available, they can be displayed using
%            ctrl+i
%          - Slightly expanded documentation
%  4.2.1   - New: Headers are stored in the 'extras' field for DICOM files
%          - Fixed: Export and Save Workspace functions remember the last path
%  4.2     - New: Added table to show standard deviations
%            New: Colormaps are now cycled through the array of functions
%            listed in prefs.colormaps
%          - New: Default radius of the circle tool is now 10% of the largest
%            image dimension
%          - New: The global struct 'svPrefs' can be used to override
%            defaultPrefs defined in this function
%          - New: Data stored in a SimpleViewer figure via
%            setappdata(hFig, 'extras'); are preserved in .svmat workspaces
%          - New: SimpleViewer workspaces now save colormaps for each subplot
%  4.1.7   - New: Toggle preview of ROI masking
%  4.1.6   - New: Support loading of all video formats supported by VideoReader
%          - New: Support custom colormaps with shift+m
%  4.1.5   - New: Grayscale animated GIFs use a more logical gray(256) colormap
%  4.1.4   - New: prefs.interpRoi can be used to interpolate ROIs only when
%            finishing an ROI for the first time (e.g. with a cubic spline) and
%            use linear interpolation otherwise (i.e. when using poly2mask)
%  4.1.3   - Fixed: Support loading of multiple animated gifs using the subplot
%            dimension
%  4.1.2   - Fixed: Improved backwards compatibility with <R2014b
%          - New: Support filenames with asterisk wildcards to load a series of
%            2D images
%          - New: 'disableVertices' preference disables the plotting of ROI
%            vertices, which can significantly improve performance (but means
%            that ROIs are not editable)
%          - Fixed: ROI plot figure does not automatically respawn after being
%            manually closed
%  4.1.1   - Fixed: resizing the figure window no longer results in black bars
%            if the image is zoomed in
%  4.1     - Fixed: figure could be displayed off-screen in certain cases
%          - New: 'CLim' can now be a 1D cell array, with CLims for each subplot
%          - Changed: XLim/YLim are saved in SimpleViewer workspaces as
%            XLim/YLim instead of XLim/YLim for consistency
%          - Changed: Workspaces are saved by default with a .svmat extension.
%            .mat files can still be saved and workspacse of either extension
%            may be loaded
%          - When loading data from images, filenames are now displayed in txt
%  4.0.9   - Fixed: crash when pasting from an empty clipboard
%          - Fixed: odd nudge behaviour if right clicking without releasing left
%            click first
%          - New: 'disableLevelAdjust' preference disables level adjustment
%            during manual window/leveling -- useful for images with a
%            meaningful zero point (e.g. positive/negative data)
%          - New: 'linkWindowLevel' preference controls which subplots (or
%            none/all) should have linked window/level adjustment
%          - New: 'blackWhiteRoi' preference
%          - Internal: InProgressRoiMove has an optimized function to improve
%            lag when drawing on a figure with a lot of points
%  4.0.8   - Bug fix: Disable ROI selection within the circle tool
%          - Bug fix: Undo while no ROIs are active no longer crashes
%          - Bug fix: Fix crash when adding a vertex in a different subplot
%          - Change: ROIs created by the circle tool no longer duplicate the
%            first point as the last point
%          - ROIs pasted from the clipboard are automatically selected
%          - Bug fix: Scroll wheel works for R2014b+
%  4.0.7   - Added support for exporting to JPEGs, PNGs, AVIs, and MP4 formats
%  4.0.6   - Fixed bug where displayed signal intensities for point ROIs are
%            offset by half a pixel
%          - SimpleViewer workspaces now save the pan/zoom state of each subplot
%  4.0.5   - New feature: Support for reading DICOM files, NIfTI files (requires
%            load_nii.m), and file type supported by imread
%  4.0.4   - New feature: Export animated GIFs using the 'e' key
%  4.0.3   - New feature: "raw data" tab shows the mean signal in each ROI in a
%            table form for quick copy/paste exporting
%  4.0.2   - New feature: "m-mode" display shows a spatial profile through time.
%            Shows up as a new tab in the ROI plot window and is enabled through
%            the pref 'showMMode'
%  4.0.1   - Compatibility fix for versions before R2013a involving gobjects
%  4.0     - Preliminary support for R2014b
%  3.0.1   - ROI tools!
%          - Logical images are now allowed
%          - Added keyboard shortcut to toggle colormap
%  2.0.2   - Double-click to use 'CLimMode' of 'auto'
%          - Must hold down "alt" to use any window/level functions
%  2.0.1   - MAJOR REVISION: Added subplot capability
%  1.0.5:  - Added mouse window/level functionality
%  1.0.4:  - Added 'strTitle' to allow figure titling
%          - Added 'scale' to set default zoom
%          - Fixed initial position so window isn't hidden behind taskbar
%            (Windows XP Default Theme)
%  1.0.3:  - Added 'space' to reset time frame
%  1.0.2:  - Viewer window can be resized and scale factor is kept (but excess
%            white space around windows is discarded)
%  1.0.1:  - Initial release
%
% Kelvin Chow (kelvin.chow@gmail.com)
% Cardiovascular MR R&D
% Siemens Healthcare
% Revision: 4.2.4  Date: 17 March 2017
%
% Copyright (c) 2017 Kelvin Chow
%
% Please do not distribute without the author's permission.

	%% --- Default settings ---
	defaultPrefs = struct(...
		'FontSize',               get(0,'DefaultAxesFontSize'), ...  % Status font size
		'homeTime',               1, ...                % Pressing spacebar goes to this time frame (can be an array for multiple spacebar presses)
		'defaultScale',           3, ...                % Initial scale factor
		'aspectRatio',            GetAspectRatio, ...   % Screen aspect ratio (used to determine automatic tiling size)
		'autoResize',             false, ...            % Keep the zoom ratio when changing time/slice?
		'roiColors',              'rgbymcwk', ...       % ROI colour order (cyclic)
		'FontName',               'Lucida Console', ... % Status font name
		'suppressRoiPlot',        true, ...             % Disable ROI signal intensity plot?
		'roiNudgeFactor',         0.125, ...            % Number of pixels moved when using shift+arrow keys
		'roiScaleFactor',         0.005, ...            % ROI scaling factor when using shrink/enlarge keys
		'staticRoi',              true, ...             % Single ROI set for all time frames in each slice?
		'autoSlicesToCine',       true, ...             % Assume a single 3D matrix is a time series instead of a multi-slice series?
		'showStaticRoiStatus',    true, ...             % Enable display of "StaticRoi is OFF" if appropriate
		'forceFieldRadius',       8, ...                % Size of nudge tool force field
		'circleRadius',           nan, ...              % Size of circle tool -- If set to NaN, this is overridden later and set to 10% of the maximum image dimension
		'circlePoints',           31, ...               % Number of points in the circle
		'enableRuler',            true, ...             % Line ROIs show their length
		'maxFrameRate',           50, ...               % Maximum frame rate at which to switch time/slice (in frames per second)
		'confirmClose',           false, ...            % Confirm window close with a dialog box
		'kbCallbacks',            struct, ...           % Additional kb callbacks (field name is key, field value is function handle). Function takes one argument -- figure handle
		'autoInterpolateRois',    false, ...            % Automatically interpolate ROIs to have n vertices
		'interpIgnoreLast',       0, ...                % If interpolating ROIs, ignore this many edges at the end
		'interpNPts',             30, ...               % Interpolate ROIs to this many vertices
		'PixelSpacing',           [1 1], ...            % Pixel spacing ([rows columns])
		'interpAlgorithm',        'linear', ...         % Interpolation between vertices within an ROI ('linear', 'spline' 'pchip')
		'transcribeRoi',          false, ...            % Draw ROIs on multiple subplots simultaneously
		'linkWindowLevel',        true, ...             % Subplots whose window/level should be linked (true, false, or cell array of subplots that should be linked)
		'showStdPlot',            true, ...             % Show +/- standard deviations on ROI data plot
		'roiLabelMode',           'none',...            % ROI label information ('none', 'number', 'area', 'signal', 'both')
		'plotShowAllSlices',      false, ...            % ???
		'disableLevelAdjust',     false, ...            % Disable level adjustment during window/level (for pos/neg data)
		'blackWhiteRoi',          false, ...            % true, false, or array of subplots with black/white ROIs
		'showMMode',              true, ...             % Show "M-Mode" image for profile through time
		'disableVertices',        false, ...            % Disable vertices on ROIs (speed boost for visualization)
		'interpRois',             'none', ...           % When to interpolate ROIs ('display', 'finish', 'always', 'none')
		'assumeMyoSaxROIs',       true, ...             % Assume short-axis myocardial ROIs when previewing ROI masks
		'colormaps',              {{@gray, @jet}}, ...  % Cell array of colormap functions to cycle through
	  'useDicomColormaps',      false);                % Use colormaps from DICOM headers, if provided via the "extras" field (e.g. via ReadDicomDir)
		
% 		'hideRoiLabels',          false, ...             % Hide ROI numeric labels (improves performance)
% 		'showRoiInfo',            false, ...            % Show ROI perimeter and area in label

	% Graphics engine is different between Mac and Windows
	if ispc
		defaultPrefs.FontSize = defaultPrefs.FontSize - 2;
	end

	% Override these default prefs if there's a global prefs struct
	global svPrefs
	
	if ~isempty(svPrefs)
		fields = fieldnames(svPrefs);
		for iField = 1:numel(fields)
			defaultPrefs.(fields{iField}) = svPrefs.(fields{iField});
		end
	end
	
	% Allow loading of files into SimpleViewer
	bIsFilename = false;
	if exist('img', 'var') && ischar(img)
		[folder, name, ext] = fileparts(img);
		if exist(img, 'file')
			if ~exist('txt', 'var') || isempty(txt)
				txt = {cat(2, name, ext)};
			end
			bIsFilename = true;
		else
			% Check if we're trying to load multiple images with a wildcard
			if any(name == '*')
				files = dir(img);
				files = files(~[files.isdir]);  % Remove folders
				if ~isempty(files)
					img = arrayfun(@(x) fullfile(folder, x.name), files, 'UniformOutput', false);
					if ~exist('txt', 'var') || isempty(txt)
						txt = {{{files.name}}};
					end
					bIsFilename = true;
				else
					error('Could not find any matching files')
				end
			end
		end
	end

% 	if exist('img', 'var') && ischar(img) && exist(img, 'file')
% 		[~, name, ext] = fileparts(img);

	if bIsFilename		
		if any(strcmpi(ext, {'.mat', '.svmat'}))
			% --- SimpleViewer .mat workspace ----------------------------------------
			ws = warning('off', 'MATLAB:load:variableNotFound');
			load(img, '-mat', 'img', 'CLim', 'roi', 'txt', 'tileSize', 'prefs', 'cmap', 'XLim', 'YLim', 'extras')
			warning(ws)
		elseif (ischar(img) && isdicom(img)) || (iscell(img) && all(cellfun(@(x) isdicom(x), img)))
			% --- DICOM files --------------------------------------------------------
			if iscell(img)
				% TODO: this will break if there's a DICOMDIR file
				info = dicominfo(img{1});
				img = cellfun(@(x) dicomread_wrapper(x), img, 'UniformOutput', false);
			else
				info = dicominfo(img);
				img = dicomread_wrapper(img);
			end

			% Pixel spacing
			% TODO: Verify that this is the same for all displayed series before continuing!
			if isfield(info, 'PixelSpacing')
				prefs.PixelSpacing = info.PixelSpacing';
			elseif isfield(info, 'SequenceOfUltrasoundRegions') && isfield(info.SequenceOfUltrasoundRegions, 'Item_1') && all(isfield(info.SequenceOfUltrasoundRegions.Item_1, {'PhysicalDeltaX', 'PhysicalDeltaY'}))
				prefs.PixelSpacing = [info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX] * 10;
			end

			% Look for RGB/truecolor data
			if isfield(info, 'ColorType') && strcmp(info.ColorType, 'truecolor') && (size(img,3) == 3)
				img = squeeze(num2cell(img, [1 2 3]));
				prefs = rmfield(prefs, 'PixelSpacing');
% 				cmap = gray(256);  % dummy
			end

			% Use custom colour map if present
			if all(isfield(info, {'RedPaletteColorLookupTableData', 'GreenPaletteColorLookupTableData', 'BluePaletteColorLookupTableData', 'RedPaletteColorLookupTableDescriptor'}))
				cmap = cat(2, info.RedPaletteColorLookupTableData, info.GreenPaletteColorLookupTableData, info.BluePaletteColorLookupTableData);
				if isrow(cmap)
					cmap = reshape(cmap, [], 3);
				end
				cmap = double(cmap) / double(max(cmap(:)));

				% Assume this is the same as Green and Blue's palette
				% The first value is the number of entries in the lookup table. The second value is the first input value mapped.
				CLim = double(info.RedPaletteColorLookupTableDescriptor(2) + [0 info.RedPaletteColorLookupTableDescriptor(1)-1]);

				% But.. MATLAB's 'painters' renderer doesn't support more than 256 colormap entries
				t = linspace(CLim(1), CLim(2), 256)';
				cmap = cat(2, interp1(CLim(1):CLim(2), cmap(:,1), t), ...
											interp1(CLim(1):CLim(2), cmap(:,2), t), ...
											interp1(CLim(1):CLim(2), cmap(:,3), t));
			else
				cmap = gray(256);
			end

			if isfield(info, 'WindowCenter') && ~isempty(info.WindowCenter) && isfield(info, 'WindowWidth') && ~isempty(info.WindowWidth)
				% Some implementation store multiple window centers; use the first one
				CLim = info.WindowCenter(1) + info.WindowWidth/2*[-1 1];
			end
			
			% Save info into workspace
			extras.headers{1}{1}{1} = info;
			
		elseif any(strcmpi(ext, {'.gz', '.nii'}))
			% --- NIfTI files --------------------------------------------------------
			if iscell(img)
				nii = cellfun(@(x) load_nii(x), img);
				img = cat(4, nii.img);
			else
				nii = load_nii(img);
				img = nii.img;
			end
			
			% If there are slices and times, then permute the slices to the subplot dimension
			if (size(img,3) > 1) && (size(img,4) > 1) && (size(img,5) == 1)
				img = permute(img, [1 2 4 5 3]);
			end

			% Make sure the header txt is the correct size if there are multiple subplots
			if (size(img,5) ~= numel(txt))
				% FIXME
% 				clear txt
				txt = repmat(txt, [1 size(img,5)]);
			end

			cmap = gray(256);
		elseif any(strcmpi(ext, {'.hdr', '.img'}))
			% --- Analyze files ------------------------------------------------------
			if ischar(img)
				files = {img};
			else
				files = img;
			end

			% Load each image (TODO: This will break if any individual image isn't a still frame)
			img = cell(1, numel(files));
			for i = 1:numel(files)
				if exist('analyze75readndim.m', 'file')
					tImg = analyze75readndim(files{i});
					if (ndims(tImg) > 4)
						tImg = squeeze(tImg);
					end
				else
					tImg = analyze75read(files{i});
				end

% 				% Check if this is complex data
% 				if strcmp(files{i}(end-7:end-4), '_MAG')
% 					if exist([files{i}(1:end-8) '_PHASE.hdr'], 'file')
% 						if exist('analyze75readndim.m', 'file')
% 							tImgPhs = analyze75readndim([files{i}(1:end-8) '_PHASE.hdr']);
% 							if (ndims(tImgPhs) > 4)
% 								tImgPhs = squeeze(tImgPhs);
% 							end
% 						else
% 							tImgPhs = analyze75read(files{i});
% 						end
% 						tImg = tImg .* exp(-(1j).*tImgPhs);
% 					end
% 				end
% 				tImg = real(tImg);

				% Strip off the path if there are multiple files
				if (numel(files) > 1)
					[~, filename, ext] = fileparts(files{i});
					filename = [filename ext];
				else
					filename = files{i};
				end
				
				if ismatrix(tImg)
					img{i} = {{tImg}};
					txt{i} = {{{filename}}};
				else
					if (ndims(tImg) > 4)
						error('>4d data not supported')
					end

					tmpCellImg = squeeze(num2cell(tImg, [1 2]));
					sz = size(tmpCellImg);

					if sz(2) == 1
						img{i} = {NestCell(tmpCellImg)};
						txt{i} = {NestCell(repmat({filename}, sz))};
					else
						img{i} = NestCell(tmpCellImg);
						txt{i} = NestCell(repmat({filename}, sz));
					end
				end
			end
			cmap = gray(256);

		elseif any(strcmpi(ext, {'.h5'}))
			% --- HDF5 files ---------------------------------------------------------
			% Each subplot shows a group within the h5 file
			if iscell(img)
				error('Loading of multiple HDF5 files not supported')
			end
			fileName = img;
			
			info       = h5info(fileName);
			
			if strcmp(info.Groups(end).Name, '/dataset')
				error('[%s] This is likely not an image HDF5 dataset', mfilename)
			end
			groupNames = {info.Groups(end).Groups.Name}';

			img        = cell(1, numel(groupNames));
			txt        = cell(1, numel(groupNames));
% 			txtComment = cell(1, numel(groupNames));
			for iGroup = 1:numel(groupNames)
				tImg  = h5read(fileName, strcat(groupNames{iGroup}, '/data'));
				tHdr  = h5read(fileName, strcat(groupNames{iGroup}, '/header'));
				tAttr = h5read(fileName, strcat(groupNames{iGroup}, '/attributes'));

				nReps = numel(unique(tHdr.repetition));
				nSets = numel(unique(tHdr.set));

				% Extract image comments
				ctTxt = cell(numel(tAttr), 1);
				for iImg = 1:numel(tAttr)
					tComments = '';
					re1 = regexp(tAttr{iImg}, '<name>GADGETRON_ImageComment</name>.*?</meta>', 'match', 'once');
					if ~isempty(re1)
						re2 = regexp(re1, '<value>(?<value>.*?)</value>', 'names');
						tComments = sprintf('%s, ', re2.value);
						tComments = [sprintf('\n'), 'Comments: ', tComments(1:end-2)];
					end
					ctTxt{iImg} = sprintf('%s\n%s%s', cat(2, name, ext), groupNames{iGroup}, tComments);
				end

				% Assume the first 2 dimensions are x,y and the remaining dimensions are
				% stored in slice/time using the nested cell format
				tImg  = squeeze(tImg);
				ctImg = permute(num2cell(tImg, [1 2]), [3 4 5 1 2]);
				if (ndims(ctTxt) > 2)
					ctTxt = permute(ctTxt, [2:ndims(tImg)-2 1]);
				end

				% Seems like there's no comment in the slices dimension??
				if (size(ctTxt,1) ~= size(ctImg,1))
					ctTxt = repmat(ctTxt, [size(ctImg,1) 1]);
				end
				
				bSortSetsReps = true;
				if bSortSetsReps
					if (nSets > 1)
						sz = size(ctImg);
						if (mod(sz(end), nSets) == 0)
							ctImg = reshape(ctImg, [sz(1:end-1) nSets sz(end)/nSets]);
							ctTxt = reshape(ctTxt, [sz(1:end-1) nSets sz(end)/nSets]);
						end
					end
				end

				img{iGroup} = NestCell(permute(ctImg, [2 1 3]));
				txt{iGroup} = NestCell(permute(ctTxt, [2 1 3]));
			end
				
% 			img = NestCell(permute(ctImg, [2 1 3]));
% 			txt = NestCell(permute(ctTxt, [2 1 3]));
% % 			img = NestCell(permute(ctImg, [3 2 1]));
% % 			txt = NestCell(permute(ctTxt, [3 2 1]));

% % % % % 				% Try to somewhat smartly reshape this
% % % % % 				if ((nReps > 1) || (nSets > 1))
% % % % % 					sz = size(tImg);
% % % % % 					if (mod(sz(3), nReps*nSets) == 0)
% % % % % 						nEco = sz(3) / (nReps*nSets);
% % % % % 						tImg = reshape(tImg, [sz([1 2]) nEco nSets nReps]);
% % % % % 						if (nEco > 1)
% % % % % 							tImg = permute(tImg, [1 2 5 4 3]);
% % % % % 						end
% % % % % 					end
% % % % % 				end
% % % % 				ctImg = squeeze(num2cell(tImg, [1 2]));
% % % 				ctImg = permute(num2cell(tImg, [1 2]), [3 4 5 1 2]);
% % % 
% % % 				bSortSetsReps = true;
% % % 				if bSortSetsReps
% % % 					if ((nReps > 1) || (nSets > 1))
% % % 						sz = size(ctImg);
% % % 						if (mod(sz(end), nSets) == 0)
% % % 							ctImg = reshape(ctImg, [sz(1:end-1) nSets sz(end)/nSets]);
% % % 							ctTxt = reshape(ctTxt, [nSets sz(end)/nSets]);
% % % 
% % % % 						if (mod(numel(ctImg), nReps*nSets) == 0)
% % % % 							ctImg = reshape(ctImg, [nSets nReps]);
% % % % 							ctTxt = reshape(ctTxt, [nSets nReps]);
% % % 						end
% % % 					end
% % % 				end
% % % 
% % % 				if (ndims(tImg) == 5)
% % % % 					img = num2cell(squeeze(num2cell(ctImg, 2)),1);
% % % 					img = num2cell(permute(num2cell(ctImg, 2), [1 3 2]),1);
% % % 					img = cellfun(@(x) x', img, 'UniformOutput', false);  % Make sure all dimensions are rows for OCD purposes
% % % 
% % % 					txt = num2cell(permute(num2cell(ctTxt, 2), [1 3 2]),1);
% % % 					txt = cellfun(@(x) x', txt, 'UniformOutput', false);  % Make sure all dimensions are rows for OCD purposes
% % % 				elseif size(ctImg,2) > 1
% % % 					img{iGroup} = cellfun(@(x) x', num2cell(ctImg,1), 'UniformOutput', false);
% % % 					txt{iGroup} = cellfun(@(x) x', num2cell(ctTxt,1), 'UniformOutput', false);
% % % 				else
% % % 					img{iGroup} = {ctImg'};
% % % 					txt{iGroup} = {ctTxt'};
% % % 				end
% % % 			end
% 			txt  = arrayfun(@(x) sprintf('%s\n%s\n%s', txt{1}, groupNames{x}, txtComment{x}), 1:numel(groupNames), 'UniformOutput', false);
% 			txt  = arrayfun(@(x) sprintf('%s\n%s\n', txt{1}, groupNames{x}), 1:numel(groupNames), 'UniformOutput', false);
% 						txt = txtComment;
% 						txt = [];
			cmap = gray(256);

% 			% Hack for when we changed the dimensions of img
% 			if numel(txt) ~= numel(img)
% % 				txt = [];
% 				nSlices = cellfun(@(x) numel(x), img);
% 				txt = arrayfun(@(x) repmat({txt(1)}, [1 x]), nSlices, 'UniformOutput', false);
% 			end
		elseif any(strcmpi(ext, {'.short', '.real', 'cplx'}))
			% --- Gadgetron Simple File Formats --------------------------------------
			if iscell(img)
				error('Loading of multiple Gadgetron files not supported')
			end
			fileName = img;

			fid = fopen(fileName);
			numdims = fread(fid, 1,       'int32'); 
			dims    = fread(fid, numdims, 'int32'); 

			switch ext
				case '.short'
					img = fread(fid, prod(dims), 'uint16'); 
				case '.real'
					img = fread(fid, prod(dims), 'float32'); 
				case '.cplx'
					img = fread(fid, 2*prod(dims), 'float32'); 
					img = complex(img(1:2:end),img(2:2:end));
			end
			fclose(fid);

			if numel(dims)>1
				img = reshape(img,dims');
			end

			cmap = gray(256);
		elseif ~isempty(imformats(ext(2:end)))
			% --- Image formats supported by imread ----------------------------------
			if ischar(img)
				files = {img};
			else
				files = img;
			end

			% Load each image (TODO: This will break if any individual image isn't a still frame)
			img = cell(1, numel(files));
			for i = 1:numel(files)
				[tImg, cmap] = imread(files{i});
				% Some formats like GIF have an embedded colormap
				if ~isempty(cmap)
					% If the colormap is all gray, convert to grayscale so that it can be window-leveled
					if all(all(diff(cmap,1,2) == 0))
						tImg = ind2gray(tImg, cmap);
						
						% See if it's a series of images (e.g. animated GIF)
						if (size(tImg,3) == 1) && (size(tImg,4) > 1)
% 							if numel(files) > 1 % fixme: why are these the same??
								img{i} = {permute(num2cell(tImg, [1 2]), [1 4 2 3])};
% 							else
% 								img{i} = {permute(num2cell(tImg, [1 2]), [1 4 2 3])};
% 							end
						else
							img{i} = tImg;
						end
						cmap = gray(256);
						CLim = [0 255];
					else
						img{i} = {permute(num2cell(tImg, [1 2]), [1 4 2 3])};
						CLim = [0 size(cmap,1)];
					end
				else%if isempty(cmap) && (size(tImg,3) == 3)
					% Some formats are read in RGB/TrueColor format
					img{i} = tImg;
% 				else
% 					error('[%s] Unsupported image format', mfilename)
				end
			end

% 			% TODO: Silly hack, fix me
% 			if (numel(files) > 1) && iscell(img{1}{1})
% 				txt = txt{1}{1};
% 			end
		else
			% Check for supported video formats
			fs = getFilterSpec(VideoReader.getFileFormats());
			vidExts = regexp(fs{1}(2:end-1), ';\*', 'split');
			vidExts = cat(2, vidExts, {'.mpeg'}); % Common variation
			
			if any(strcmp(ext, vidExts))
				% --- Video formats supported by VideoReader -----------------------------
				v = VideoReader(img);
				cImg = {};
				while hasFrame(v)
					cImg{end+1} = readFrame(v);
				end
				close(v);
				img = cImg;

				% QuickTime videos have an odd skew issue
				if strcmp(ext, '.mov')
					sz = size(img{end});

					for iFrame = 1:numel(img)
						for iRow = 1:sz(1)
							ind = (1:sz(1))-iRow-1;
							ind = mod(ind-1, sz(1)) + 1;
							img{iFrame}(iRow,:,:) = img{iFrame}(iRow,ind,:);
						end
					end
					
					% Disable support for RGB because it's not known how to support it now
					img = cat(4, img{:});
					img = img(:,:,1,:);
					disp(sprintf('Stored frame rate is: %1.1f fps', v.FrameRate))
				else
					% If RGB channels are the same for all time, then it's grayscale and we can just use one channel
					if all(cellfun(@(x) all(~reshape(diff(x,1,3), 1, [])), img))
						img = cat(4, img{:});
						img = img(:,:,1,:);
					end
				end
				cmap = gray(256);
			else
				error('[%s] Unsupported file type', mfilename)
			end
		end
	else
		if ~exist('cmap', 'var') || isempty(cmap)
			cmap = gray(256);
		end
	end
	
	% Display some info if no arguments are passed
	if (nargin == 0) || isempty(img)
		disp(sprintf('Syntax: %s(img, CLim, roi, txt, tileSize, prefs)', mfilename))

		disp('Preference fields:')
		f = fieldnames(defaultPrefs);
		nCols = 5;
		nRows = ceil(numel(f)/nCols);

		if numel(f) < nRows*nCols
			f{nRows*nCols} = [];
		end
		
		if exist('CreateTextTable.m', 'file')
			disp(CreateTextTable(reshape(cat(1, cell(nRows,1), f), nRows, nCols+1)))
		else
			disp(reshape(f, nRows, nCols))
		end
		return
	end

	% Allow the prefs struct to be any argument
	if exist('CLim', 'var') && isstruct(CLim)
		prefs = CLim;
		clear CLim
	end

	if exist('roi', 'var') && isstruct(roi)
		prefs = roi;
		clear roi
	end
	
	if exist('txt', 'var') && isstruct(txt)
		prefs = txt;
		clear txt
	end

	if exist('tileSize', 'var') && isstruct(tileSize)
		prefs = tileSize;
		clear tileSize
	end

	if exist('XLim', 'var') && isstruct(XLim)
		prefs = XLim;
		clear XLim
	end

	if exist('YLim', 'var') && isstruct(YLim)
		prefs = YLim;
		clear YLim
	end

	if exist('cmap', 'var') && isstruct(cmap)
		prefs = cmap;
		clear cmap
	end

	if ~exist('prefs', 'var')
		prefs = struct;
	end

	%% 1. Validate inputs
	% --- Image data -------------------------------------------------------------
	if issparse(img)
		img = full(img);
	end
	
	if islogical(img)
		img = uint8(img);
	end

	if isnumeric(img)
		if ndims(img) < 2 || ndims(img) > 5
			error('"img" should be 3, 4, or 5 dimensional when an array')
		end
	elseif iscell(img)
		if (numel(img) ~= length(img)) || isempty(img)
			error('"img" should a nested 1 dimensional cell arrays')
		end
	else
		error('"img" should be an array, a cell array, or a file name')
	end

	if isnumeric(img)
		imgData = cell(1, size(img,5));
		for iPlot = 1:size(img,5)
			imgData{iPlot} = cell(1, size(img,3));
			for iSlice = 1:size(img,3)
				imgData{iPlot}{iSlice} = cell(1, size(img,4));
				for iTime = 1:size(img,4)
					imgData{iPlot}{iSlice}{iTime} = img(:,:, iSlice, iTime, iPlot);
				end
			end
		end
	elseif iscell(img)
		% 1D cell array of images
		if all(cellfun(@(x) isnumeric(x), img))
			if (isfield(prefs, 'autoSlicesToCine') && prefs.autoSlicesToCine) || defaultPrefs.autoSlicesToCine
				imgData = {{img}};
			else
				imgData = {cellfun(@(x) {x}, img, 'UniformOutput', false)};
			end
		else
			imgData = img;
		end
	end

	% Check for complex data
	bWarnComplex = false;
	for iPlot = 1:numel(imgData)
		for iSlice = 1:numel(imgData{iPlot})
			for iTime = 1:numel(imgData{iPlot}{iSlice})
				if ~isreal(imgData{iPlot}{iSlice}{iTime})
					imgData{iPlot}{iSlice}{iTime} = abs(imgData{iPlot}{iSlice}{iTime});
					bWarnComplex = true;
				end
			end
		end
	end
	if bWarnComplex
		disp('Warning: Converted complex img to abs')
	end

	% Allow txt to be a 1D matrix for each subplot
	if exist('txt', 'var') && iscell(txt) && all(cellfun(@(x) ischar(x), txt))
		txt = cellfun(@(x) {{x}}, txt, 'UniformOutput', false);
	end

	% Allow txt to be a 1D matrix for each subplot
% 	if exist('txt', 'var') && iscell(txt) && all(cellfun(@(x) ischar(x), txt))
% 		txt = cellfun(@(x) {{x}}, txt, 'UniformOutput', false);
% 	end

	% --- ROI, text label data ---------------------------------------------------
	bDummyRoi = false;
	if ~exist('roi', 'var') || isempty(roi)
		bDummyRoi = true;
		roiData = cell(1, numel(imgData));
		for iPlot = 1:numel(imgData)
			roiData{iPlot} = cell(1, numel(imgData{iPlot}));
			for iSlice = 1:numel(imgData{iPlot})
				roiData{iPlot}{iSlice} = cell(1, numel(imgData{iPlot}{iSlice}));
				for iTime = 1:numel(imgData{iPlot}{iSlice})
					roiData{iPlot}{iSlice}{iTime} = cell(1,0);
				end
			end
		end
	else
		roiData = roi;
	end

	bDummyText = false;
	if ~exist('txt', 'var') || isempty(txt)
		bDummyText = true;
		txtData = cell(1, numel(imgData));
		for iPlot = 1:numel(imgData)
			txtData{iPlot} = cell(1, numel(imgData{iPlot}));
			for iSlice = 1:numel(imgData{iPlot})
				txtData{iPlot}{iSlice} = cell(1, numel(imgData{iPlot}{iSlice}));
				for iTime = 1:numel(imgData{iPlot}{iSlice})
					txtData{iPlot}{iSlice}{iTime} = cell(1,0);
				end
			end
		end
	else
		if iscell(txt) && all(cellfun(@(x) ischar(x), txt))
			if (isfield(prefs, 'autoSlicesToCine') && prefs.autoSlicesToCine) || defaultPrefs.autoSlicesToCine
				txtData = {{txt}};
			else
				txtData = {cellfun(@(x) {x}, txt, 'UniformOutput', false)};
			end
		else
			txtData = txt;
		end
	end

	%% 1.5 Single-slice cines
	% Note: do this early, and only for images (or dummy txt/roi data) PRIOR to the dimensions check
	if (isfield(prefs, 'autoSlicesToCine') && prefs.autoSlicesToCine) || defaultPrefs.autoSlicesToCine
		% If there is only a single series of images, assume this is a time-series
		% (movie) and transpose this to the time dimension
		bSingleSliceCine = true;
		for iPlot = 1:numel(imgData)
			if ~isempty(imgData{iPlot}) && any(cellfun(@(x) numel(x) ~= 1, imgData{iPlot}))
				bSingleSliceCine = false;
			end
		end

		if bSingleSliceCine
			oldImgData = imgData;  imgData = cell(1, numel(oldImgData));
			for iPlot = 1:numel(oldImgData)
				imgData{iPlot} = {cell(1, numel(oldImgData{iPlot}))};
				for iTime = 1:numel(oldImgData{iPlot})
					imgData{iPlot}{1}{iTime} = oldImgData{iPlot}{iTime}{1};
				end
			end

			% Convert only if it's dummy text data
			if bDummyText || numel(oldImgData{1}) == numel(txtData{1})
				oldTxtData = txtData;  txtData = cell(1, numel(oldTxtData));
				for iPlot = 1:numel(oldTxtData)
					txtData{iPlot} = {cell(1, numel(oldTxtData{iPlot}))};
					for iTime = 1:numel(oldTxtData{iPlot})
						txtData{iPlot}{1}{iTime} = oldTxtData{iPlot}{iTime}{1};
					end
				end
			end

			% Convert only if it's dummy ROI data or if ROIs were saved as single slice too
			if bDummyRoi || numel(oldImgData{1}) == numel(roiData{1})
				oldRoiData = roiData;  roiData = cell(1, numel(oldImgData));
				for iPlot = 1:numel(oldRoiData)
					roiData{iPlot} = {cell(1, numel(oldRoiData{iPlot}))};
					for iTime = 1:numel(oldRoiData{iPlot})
						roiData{iPlot}{1}{iTime} = oldRoiData{iPlot}{iTime}{1};
					end
				end
			end
		end
	end
	
	% Check that roi and txt have the same dimensions as img
	if ~iscell(roiData)
		error('''roi'' should be a nested cell');
	elseif ~iscell(txtData)
		error('''txt'' should be a nested cell');
	end
	if numel(roiData) ~= numel(imgData)
		error('''roi'' should be a nested cell with the same number of subplots/slices/times as ''img''');
	elseif numel(txtData) ~= numel(imgData)
% 		error('''txt'' should be a nested cell with the same number of subplots/slices/times as ''img''');
	end
	for iPlot = 1:numel(imgData)
		if ~iscell(roiData{iPlot})
			error('''roi'' should be a nested cell with the same number of subplots/slices/times as ''img''');
		elseif ~iscell(txtData{iPlot})
			error('''txt'' should be a nested cell with the same number of subplots/slices/times as ''img''');
		end
		if numel(roiData{iPlot}) ~= numel(imgData{iPlot})
			error('''roi'' should be a nested cell with the same number of subplots/slices/times as ''img''');
		elseif numel(txtData{iPlot}) ~= numel(imgData{iPlot})
% 			error('''txt'' should be a nested cell with the same number of subplots/slices/times as ''img''');
		end
		for iSlice = 1:numel(imgData{iPlot})
			if ~iscell(roiData{iPlot}{iSlice})
				error('''roi'' should be a nested cell with the same number of subplots/slices/times as ''img''');
			elseif ~iscell(txtData{iPlot}{iSlice})
% 				error('''txt'' should be a nested cell with the same number of subplots/slices/times as ''img''');
			end
			if numel(roiData{iPlot}{iSlice}) ~= numel(imgData{iPlot}{iSlice})
% 				error('Not this time, Mr. Bond...');
				error('''roi'' should be a nested cell with the same number of subplots/slices/times as ''img''');
			elseif numel(txtData{iPlot}{iSlice}) ~= numel(imgData{iPlot}{iSlice})
				% This is not an error, as we're allowing a single text label for the all times in a slice
% 				error('''txt'' should be a nested cell with the same number of subplots/slices/times as ''img''');
			end
		end
	end

	% Allow [] instead of {1x0 cell} for a given time frame
	for iPlot = 1:numel(roiData)
		for iSlice = 1:numel(roiData{iPlot})
			for iTime = 1:numel(roiData{iPlot}{iSlice})
				if isempty(roiData{iPlot}{iSlice}{iTime}) && isnumeric(roiData{iPlot}{iSlice}{iTime})
					roiData{iPlot}{iSlice}{iTime} = cell([1 0]);
				end
			end
		end
	end

%% 1.5 Single-slice cines
% % % % NOTE: Old behaviour.  'autoSlicesToCine' now converts ONLY the image data.
% % % % Uncomment if necessary to revert:
% % % 	if prefs.autoSlicesToCine
% % % % 		% If there is only a single series of images, assume this is a time-series
% % % % 		% (movie) and transpose this to the time dimension
% % % % 		bSingleSliceCine = true;
% % % % 		for iPlot = 1:numel(imgData)
% % % % 			if any(cellfun(@(x) numel(x) ~= 1, imgData{iPlot}))
% % % % 				bSingleSliceCine = false;
% % % % 			end
% % % % 		end
% % % 
% % % 		if bSingleSliceCine
% % % % 			oldImgData = imgData;  imgData = cell(1, numel(oldImgData));
% % % 			oldTxtData = txtData;  txtData = cell(1, numel(oldTxtData));
% % % 			oldRoiData = roiData;  roiData = cell(1, numel(oldImgData));
% % % 			for iPlot = 1:numel(oldTxtData)
% % % % 				imgData{iPlot} = {cell(1, numel(oldImgData{iPlot}))};
% % % 				txtData{iPlot} = {cell(1, numel(oldTxtData{iPlot}))};
% % % 				roiData{iPlot} = {cell(1, numel(oldRoiData{iPlot}))};
% % % 				for iTime = 1:numel(oldTxtData{iPlot})
% % % % 					imgData{iPlot}{1}{iTime} = oldImgData{iPlot}{iTime}{1};
% % % 					txtData{iPlot}{1}{iTime} = oldTxtData{iPlot}{iTime}{1};
% % % 					roiData{iPlot}{1}{iTime} = oldRoiData{iPlot}{iTime}{1};
% % % 				end
% % % 			end
% % % 		end
% % % 	end

	%% 4. Create figure
	sz = size(imgData{1}{1}{1}); % May need to tweak this to be robust to the first image not existing
	if isnan(defaultPrefs.circleRadius)
		defaultPrefs.circleRadius = round(max(sz)/10);
	end
	prefs = GeneratePrefs(defaultPrefs, prefs);
	prefs.roiColors = ConvertColorsToRGB(prefs.roiColors);

	% Sanity check to avoid plotting way too many points
	for iPlot = 1:numel(roiData)
		for iSlice = 1:numel(roiData{iPlot})
			for iTime = 1:numel(roiData{iPlot}{iSlice})
				if ~isempty(roiData{iPlot}{iSlice}{iTime}) && iscell(roiData{iPlot}{iSlice}{iTime})
					if any(cellfun(@(x) numel(x), roi{iPlot}{iSlice}{iTime}) > 200)
						prefs.disableVertices = true;
					end
				end
			end
		end
	end
	
	if ~exist('CLim', 'var') || ~((isnumeric(CLim) && (numel(CLim) == 2)) || (iscell(CLim) && ~isempty(CLim)))
		CLim = [];
	end

	if exist('tileSize', 'var') && isnumeric(tileSize) && (numel(tileSize) == 2)
		if numel(imgData) > prod(tileSize)
			disp('tileSize is too small for the number of subplots.  Overriding...')
			clear tileSize
		end
	end
	
	if ~exist('tileSize', 'var') || ~isnumeric(tileSize) || numel(tileSize) ~= 2
		rows     = round(sqrt(numel(imgData)/prefs.aspectRatio));
		cols     = ceil(numel(imgData)/rows);
		tileSize = [rows cols];
	end

	nRows = tileSize(1);  height = 1/nRows;
	nCols = tileSize(2);  width  = 1/nCols;

	hFig  = gcf;
	clf(hFig, 'reset');

	% Status text
	setappdata(hFig, 'staticRoi',         prefs.staticRoi);
	setappdata(hFig, 'strMode',           'Normal');
	setappdata(hFig, 'strError',          '');

	setappdata(hFig, 'prefs',             prefs);
	setappdata(hFig, 'roiPlot',           struct('hFig', [], 'hAxes', []));
	setappdata(hFig, 'clipboardRoi',      []);
% 	setappdata(hFig, 'selectedSubplot', 1);
	setappdata(hFig, 'selectedRoi',       []);
	setappdata(hFig, 'hideRois',          false);
	setappdata(hFig, 'fcnUpdateRoiData',  @UpdateRoiData);
	setappdata(hFig, 'inProgressRoiAxis', []);

	% Function handles
	fcnHandles = struct('SetTimeSlice',          @SetTimeSlice, ...
	                    'SetTimeSliceImageOnly', @SetTimeSliceImageOnly, ...
	                    'UpdateRoiData',         @UpdateRoiData, ...
	                    'UpdateRoiHandles',      @UpdateRoiHandles, ...
	                    'UpdateAllRoiHandles',   @UpdateAllRoiHandles, ...
	                    'SetNormalMode',         @SetNormalMode, ...
	                    'ToggleWindowLevelMode', @ToggleWindowLevelMode, ...
	                    'ToggleRoiDrawMode',     @ToggleRoiDrawMode, ...
	                    'SetFigureScale',        @SetFigureScale, ...
	                    'MoveCircle',            @MoveCircle, ...
	                    'UpdateStatusText',      @UpdateStatusText, ...
	                    'MouseDownSelect',       @MouseDownSelect, ...
	                    'SelectSubplot',         @SelectSubplot, ...
	                    'StartRoiMove',          @StartRoiMove, ...
	                    'VertexButtonDown',      @VertexButtonDown, ...
	                    'AddVertexFromEdge',     @AddVertexFromEdge, ...
	                    'FreezeColormap',        @FreezeColormap, ...
	                    'UnFreezeColormap',      @UnFreezeColormap, ...
	                    'AutoWindowLevel',       @AutoWindowLevel, ...
	                    'ToggleCircle',          @ToggleCircle, ...
	                    'HideExtras',            @HideExtras, ...
	                    'UpdateMMode',           @UpdateMMode, ...
	                    'UpdateRawTable',        @UpdateRawTable, ...
	                    'AdjustZoom',            @AdjustZoom);

	setappdata(hFig, 'fcnHandles', fcnHandles)

	% Keypress/mouse handlers
	set(hFig, 'KeyPressFcn',         @CatchKeyPress);
	set(hFig, 'CloseRequestFcn',     @CloseSimpleViewer);

	set(hFig, 'WindowButtonDownFcn',  @MouseDownSelect);
	set(hFig, 'WindowButtonUpFcn',    '');
	set(hFig, 'WindowScrollWheelFcn', @MouseScroll);
	set(hFig, 'ResizeFcn',            @(h, e)AdjustZoom(h, e, 1));

	% hAxes has rows/columns REVERSED that of the subplot so that 1D indexing works
	hAxes = zeros(nCols, nRows);
	for iPlot = 1:numel(imgData)
		% Axes
		[iCol, iRow] = ind2sub([nCols nRows], iPlot);
		left   = (iCol - 1)     * width;
		bottom = (nRows - iRow) * height;
		CurrentAxes = subplot('Position', [left bottom width height]);
		hAxes(iPlot) = CurrentAxes;

		if isempty(imgData{iPlot})
			axis off
			continue
		end

		% Images
		try
			hImage = imagesc(imgData{iPlot}{1}{1});
		catch
			hImage = imagesc([]);
		end
		if exist('CLim', 'var') && ~isempty(CLim)
			if iscell(CLim)
				if ~isempty(CLim{iPlot})
					set(CurrentAxes, 'CLim', CLim{iPlot});
				end
			else
				set(CurrentAxes, 'CLim', CLim);
			end
		end
		
		if iscell(cmap)
			colormap(CurrentAxes, cmap{iPlot});
		end
		% Performance issues on R2016b
% 		axis image off
		set(CurrentAxes, 'DataAspectRatio', [1 1 1])
		set(CurrentAxes, 'Visible',         'off')
		hold on

		if (prefs.PixelSpacing(2) / prefs.PixelSpacing(1)) ~= 1
			axis normal
			set(CurrentAxes, 'DataAspectRatio', [prefs.PixelSpacing(:)' 1])
		end

		% Don't allow automatic resizing of limits.  Improves performance but is
		% also necessary in R2014b+ as otherwise moving the circle tool to the edge
		% will automatically resize the axis to always show all points
		set(CurrentAxes, 'XLimMode', 'manual')
		set(CurrentAxes, 'YLimMode', 'manual')
		set(CurrentAxes, 'ZLimMode', 'manual')

		% Set zoom (xlim/ylim) if specified
		if exist('XLim', 'var') && ~isempty(XLim)
			if iscell(XLim)
				set(CurrentAxes, 'XLim', XLim{iPlot})
			else
				set(CurrentAxes, 'XLim', XLim)
			end
		end

		if exist('YLim', 'var') && ~isempty(YLim)
			if iscell(YLim)
				set(CurrentAxes, 'YLim', YLim{iPlot})
			else
				set(CurrentAxes, 'YLim', YLim)
			end
		end

		% ROIs
		% The only thing we need to add is 'roi'.  UpdateRoiData will take care of the rest		
		roiHandles = [];

		% TODO: This shouldn't be necessary... UpdateRoiData should handle this case itself!
		emptyCell = cell(1, numel(imgData{iPlot}));
		for iSli = 1:numel(imgData{iPlot})
			emptyCell{iSli} = cell(1, numel(imgData{iPlot}{iSli}));
			for iTime = 1:numel(imgData{iPlot}{iSli})
				emptyCell{iSli}{iTime} = cell(1,0);
			end
		end

		setappdata(CurrentAxes, 'roi',           roiData{iPlot});
		setappdata(CurrentAxes, 'roiMask',       emptyCell);
		setappdata(CurrentAxes, 'roiData',       emptyCell);
		setappdata(CurrentAxes, 'roiStd',        emptyCell);
		setappdata(CurrentAxes, 'roiPerimeter',  emptyCell);
		setappdata(CurrentAxes, 'roiArea',       emptyCell);
		setappdata(CurrentAxes, 'roiProfile',    emptyCell);
		setappdata(CurrentAxes, 'roiHandles',    roiHandles);
		setappdata(CurrentAxes, 'inProgressRoi', []);

		% Store data
		setappdata(CurrentAxes, 'hImage',       hImage);
% % % 		setappdata(CurrentAxes, 'hText',        hText);
		setappdata(CurrentAxes, 'imgData',      imgData{iPlot});
		setappdata(CurrentAxes, 'txtData',      txtData{iPlot});
		setappdata(CurrentAxes, 'currentTime',  1);
		setappdata(CurrentAxes, 'currentSlice', 1);

		if verLessThan('matlab','8.4.0')
			set(CurrentAxes, 'Drawmode', 'fast');
		else
			set(CurrentAxes, 'SortMethod', 'childorder');
		end
	end
	setappdata(hFig, 'hAxes',           hAxes);
	setappdata(hFig, 'lastCurrentAxes', hAxes(iPlot));
	
	% Should this be done for all subplots?
% 		UpdateRoiData(hFig, iRoi);

	% This is a terribly slow way of doing this!
	for iPlot = 1:numel(imgData)
		for iSli = 1:numel(roiData{iPlot})
			setappdata(CurrentAxes, 'currentTime', iSli);
			for iTime = 1:numel(roiData{iPlot}{iSli})
				for iRoi = 1:numel(roiData{iPlot}{iSli}{iTime})
					UpdateRoiData(hAxes(iPlot), iRoi, iTime);
				end
			end
		end
	end
	setappdata(CurrentAxes, 'currentTime',  1);

	for iRoi = 1:numel(roiData{iPlot}{1}{1})
		for iPlot = 1:numel(imgData)
			UpdateRoiData(hAxes(iPlot), iRoi);
		end
	end
	
	% An invisible axes hosts the "axes selection (indicator) box" and the status text
	hInvisAxis = axes('Position', [0 0 1 1]);
	set(hInvisAxis, 'HitTest', 'off')
	if ~verLessThan('matlab','8.4.0')
		set(hInvisAxis, 'PickableParts', 'none')
	end
	hold on

	
	% Text
	for iPlot = 1:numel(imgData)
		[iCol, iRow] = ind2sub([nCols nRows], iPlot);
		left   = (iCol - 1)     * width;
		bottom = (nRows - iRow) * height;

		hText = text(left+0.01,bottom+height-0.005, '', ...
		             'Color',             'w', ...
	               'BackgroundColor',     'k', ...
		             'VerticalAlignment', 'top', ...
		             'FontName',          prefs.FontName, ...
		             'FontSize',          prefs.FontSize, ...
		             'Interpreter',       'none');

		if ~verLessThan('matlab','8.4.0')
			hText.Margin        = 2;
			hText.FontSmoothing = 'off';
		end
		if ~isempty(imgData{iPlot})
			if ischar(txtData{iPlot}{1})
				strText = ComposeText(txtData{iPlot}{1}, 1, length(imgData{iPlot}), 1, length(imgData{iPlot}{1}));
			else
				strText = ComposeText(txtData{iPlot}{1}{1}, 1, length(imgData{iPlot}), 1, length(imgData{iPlot}{1}));
			end
		else
			strText = '';
		end
		set(hText, 'String', strText);

		setappdata(hAxes(iPlot), 'hText',        hText);
	end

	hStatusText = text(0.01, 0.005, '', ...
	                   'Color',               'w', ...
	                   'BackgroundColor',     'k', ...
	                   'HorizontalAlignment', 'left', ...
	                   'VerticalAlignment',   'bottom', ...
	                   'FontName',            prefs.FontName, ...
	                   'FontSize',            prefs.FontSize, ...
	                   'Interpreter',         'none');
		if ~verLessThan('matlab','8.4.0')
			hStatusText.Margin        = 2;
			hStatusText.FontSmoothing = 'off';
		end
	axesPosition = get(hAxes(1), 'Position');
	if numel(hAxes) > 1
		hAxesBox = plot(axesPosition(1) + [0 axesPosition(3) axesPosition(3) 0 0] + 0.001*[1 -1 -1 1 1], ...
		                axesPosition(2) + [0 0 axesPosition(4) axesPosition(4) 0] + 0.001*[1 1 -1 -1 1], '-', 'Color', 'r', 'LineWidth', 1.5);
	else
		hAxesBox = [];
	end
	xlim([0 1]), ylim([0 1]);
	set(hInvisAxis, 'Visible', 'off');
	setappdata(hFig, 'hAxesBox',        hAxesBox);
	setappdata(hFig, 'hStatusText',     hStatusText);
	UpdateStatusText(hFig);

	if ~iscell(cmap)
		colormap(cmap)
	end

	% Figure tweaking
	set(hFig, 'Toolbar',     'none');
	set(hFig, 'MenuBar',     'none');
	set(hFig, 'Color',       'k');

	set(hFig, 'CurrentAxes', hAxes(1));  % Size according to the first image
	
	minSize = 2; % Minimum side length [inches]
	SetFigureScale(hFig, prefs.defaultScale, minSize);
% 	set(hFig, 'CurrentAxes', hAxes(find(hAxes,1,'last')));  % Go to last axes

	UpdateAllRoiHandles(hFig);
% % 	% Cheating again... :/
% % 	for i = 1:10
% % 		UpdateRoiData(gcf, i);
% % 	end

	setappdata(hFig, 'hSV3D', []);
	setappdata(hFig, 'hSV2D', []);

	CreateCustomPointers(hFig)
	% To fix: "text with background color (including data tips) and text displayed on image, patch, or surface objects is not visible when using OpenGL renderer."
	opengl('OpenGLBitmapZbufferBug',1)

	% Z-buffer seems faster, with no observable side-effects... so far.
	if verLessThan('matlab','8.4.0')
		set(hFig, 'Renderer', 'zbuffer')
	else
		set(hFig, 'Renderer', 'opengl')
	end
	
	% Re-load extras
	if exist('extras', 'var') && ~isempty(extras)
		setappdata(hFig, 'extras', extras);
	end
	
	% Set up menu items
	% hFig = figure; imagesc(ones(100))
% 	uiSV          = uimenu(hFig,    'Label', 'SimpleViewer');
% 	uiModes       = uimenu(uiSV,    'Label', 'Modes');
% 	uiWindowLevel = uimenu(uiModes, 'Label', 'Window/Level', 'Callback', @ToggleWindowLevelMode);
% 	uiPanZoom     = uimenu(uiModes, 'Label', 'Pan/Zoom');
% 	uiDrawRoi     = uimenu(uiModes, 'Label', 'Draw ROI');

end
	%% ------------------------------------------------------------
	function CatchKeyPress(hFig, evt)
	% Dispatch key presses to the appropriate functions
		hAxis = get(hFig, 'CurrentAxes');
		prefs = getappdata(hFig, 'prefs');

		% Allow prefs.kbCallbacks to override existing functions
		if isfield(prefs.kbCallbacks, evt.Key)
			fcnHandle = prefs.kbCallbacks.(evt.Key);
			if ~isempty(fcnHandle)
				fcnHandle(hFig);
			end
			return
		end

		if any(cellfun(@(x) strcmp(x, 'control'), evt.Modifier))
			switch(evt.Key)
				case 'space'
					StartStopMovie(hFig);
				case 'i'
					CreateInfoFigure(hFig);
			end
			return
		end

		if any(cellfun(@(x) strcmp(x, 'shift'), evt.Modifier))
			switch(evt.Key)
				case 'leftarrow'
					KbMoveRoi(hFig, -prefs.roiNudgeFactor, 0)
				case 'rightarrow'
					KbMoveRoi(hFig, prefs.roiNudgeFactor, 0)
				case 'leftbracket'
					ScaleCircle(hFig, 1-prefs.roiScaleFactor*5);
				case 'rightbracket'
					ScaleCircle(hFig, 1+prefs.roiScaleFactor*5);
				case 'uparrow'
					KbMoveRoi(hFig, 0, -prefs.roiNudgeFactor);
				case 'downarrow'
					KbMoveRoi(hFig, 0,  prefs.roiNudgeFactor);
				case 'comma'
					KbScaleRoi(hFig, 1-prefs.roiScaleFactor*5);
				case 'period'
					KbScaleRoi(hFig, 1+prefs.roiScaleFactor*5);
				case {'delete', 'backspace'}
					DeleteRois(hFig, []);
				case 'c'
					CopyRoi(hFig, true);
				case 'v'
					PasteRoi(hFig, true);
				case 'o'
					AutoMakeCircle(hFig);
				case 's'
					SaveSvWorkspace(hFig);
				case 'i'
					InterpolateImages(hFig);
				case 'm'
					ToggleColormap(hFig, true); % Custom colormap using text entry
					figure(hFig); % Bring focus back after command window entry
				case 'y'
					% Testing
					ToggleGridfitRoi(hFig);
				case 'backslash'
					SetFigureScale(hFig, 1)
				case 'space'
					StartStopMovie(hFig);
					return
			end
			% Shift-tab switches subplots
			if ~strcmp(evt.Key, 'tab')
				return
			end
		end

		if any(cellfun(@(x) strcmp(x, 'alt'), evt.Modifier))
			% Alt + Number keys (1 through 9) select subplots
			if length(evt.Key) == 1 && (evt.Key >= 49) && (evt.Key <= 57)
				hAxes = getappdata(hFig, 'hAxes');

				iSubplot = uint8(evt.Key - 48);
				if iSubplot <= numel(hAxes)
					set(hFig, 'CurrentAxes', hAxes(iSubplot))
					SelectSubplot(hFig, hAxes(iSubplot))
				end
			end
			
			if strcmp(evt.Key, 'c')
				currRenderer = get(hFig, 'Renderer');
				set(hFig, 'Renderer', 'painters')
				set(hFig, 'PaperPositionMode', 'auto')
				drawnow
				print -dbitmap -clipboard
				set(hFig, 'Renderer', currRenderer)
				disp('Figure copied to clipboard.')
			end
		else
			% Number keys (1 through 9) select their corresponding ROI
			if length(evt.Key) == 1 && (evt.Key >= 49) && (evt.Key <= 57)
				roiHandles     = getappdata(hAxis, 'roiHandles');
				oldSelectedRoi = getappdata(hFig,  'selectedRoi');
				iRoi = uint8(evt.Key - 48);
				if (iRoi <= numel(roiHandles)) && ~isempty(roiHandles(iRoi).center) && (roiHandles(iRoi).center ~= 0)
					if (iRoi == oldSelectedRoi)
						% Unselect ROI
						set(roiHandles(oldSelectedRoi).line,     'LineStyle', '-');
						set(roiHandles(oldSelectedRoi).vertices, 'Marker',    '.');
						if ismac
							set(roiHandles(oldSelectedRoi).vertices, 'MarkerSize', 12)
						end
						setappdata(hFig, 'selectedRoi', []);
					else
						SelectRoi(hFig, hAxis, hAxis, roiHandles(iRoi).center)
					end
				end
			elseif length(evt.Key) == 1 && (evt.Key == 48)
				% Select all ROIs
				roiHandles = getappdata(hAxis, 'roiHandles');
				if ~isempty(roiHandles)
					roi        = getappdata(hAxis, 'roi');
					nonEmptyRois = cellfun(@(x) ~isempty(x), roi{getappdata(hAxis, 'currentSlice')}{getappdata(hAxis, 'currentTime')});
					SelectRoi(hFig, hAxis, hAxis, [roiHandles(nonEmptyRois).center])
				end
			end
		end

		switch(evt.Key)
			case 'leftarrow'
				SetTimeSlice(hFig, getappdata(hAxis, 'currentTime')-1, getappdata(hAxis, 'currentSlice'), false, any(cellfun(@(x) strcmp(x, 'alt'), evt.Modifier)));
			case 'rightarrow'
% 				SetTimeSlice(hFig, getappdata(hAxis, 'currentTime')+1, getappdata(hAxis, 'currentSlice'));
				SetTimeSlice(hFig, getappdata(hAxis, 'currentTime')+1, getappdata(hAxis, 'currentSlice'), false, any(cellfun(@(x) strcmp(x, 'alt'), evt.Modifier)));
% 				SetTimeSliceStruct(hFig, getappdata(hAxis, 'currentTime')+1, getappdata(hAxis, 'currentSlice'));
			case 'uparrow'
				SetTimeSlice(hFig, getappdata(hAxis, 'currentTime'),   getappdata(hAxis, 'currentSlice')-1, false, any(cellfun(@(x) strcmp(x, 'alt'), evt.Modifier)));
			case 'downarrow'
				SetTimeSlice(hFig, getappdata(hAxis, 'currentTime'),   getappdata(hAxis, 'currentSlice')+1, false, any(cellfun(@(x) strcmp(x, 'alt'), evt.Modifier)));
			case {'numpad0', 'slash'}
				lastTime  = getappdata(hAxis, 'lastTime');
				lastSlice = getappdata(hAxis, 'lastSlice');
				if ~isempty(lastTime) && ~isempty(lastSlice)
					SetTimeSlice(hFig, lastTime, lastSlice);
				end
			case 'space'
				GotoHomeTime(hFig);
			case 'tab'
				hAxes = getappdata(hFig, 'hAxes');
				CurrentAxes = get(hFig, 'CurrentAxes');
				CurrentSubplot = find(hAxes == CurrentAxes);
				numSubplots = sum(hAxes(:) ~= 0);

				if any(cellfun(@(x) strcmp(x, 'shift'), evt.Modifier))
					iSubplot = mod(CurrentSubplot-1-1,numSubplots)+1;
				else
					iSubplot = mod(CurrentSubplot+1-1,numSubplots)+1;
				end

				set(hFig, 'CurrentAxes', hAxes(iSubplot))
				SelectSubplot(hFig, hAxes(iSubplot))
			case 'escape'
				CloseSimpleViewer(hFig)
			case {'add', 'equal'}
				SetFigureScale(hFig, CalculateCurrentScale(hFig)*1.25)
			case {'subtract', 'hyphen'}
				SetFigureScale(hFig, CalculateCurrentScale(hFig)/1.25)
			case 'backslash'
				SetFigureScale(hFig, prefs.defaultScale)
			case 'backquote'
				% Select all ROIs
				roiHandles = getappdata(hAxis, 'roiHandles');
				if ~isempty(roiHandles)
					roi        = getappdata(hAxis, 'roi');
					nonEmptyRois = cellfun(@(x) ~isempty(x), roi{getappdata(hAxis, 'currentSlice')}{getappdata(hAxis, 'currentTime')});
					SelectRoi(hFig, hAxis, hAxis, [roiHandles(nonEmptyRois).center])
				end
			case {'delete', 'backspace'}
				if ~isempty(getappdata(hAxis, 'inProgressRoi'))
					InProgressRoiDeleteVertex(hFig, evt)
				else
					if ~isempty(getappdata(hFig, 'selectedRoi'))
						DeleteRois(hFig, getappdata(hFig, 'selectedRoi'))
					end
				end
			case 'leftbracket'
				ScaleCircle(hFig, 1-prefs.roiScaleFactor);
			case 'rightbracket'
				ScaleCircle(hFig, 1+prefs.roiScaleFactor);
			case 'e'
				ExportImage(hFig);
			case 'm'
				ToggleColormap(hFig);
			case 'd'	
				ToggleRoiDrawMode(hFig);
			case 'w'
				ToggleWindowLevelMode(hFig);
			case 's'
				if ~any(cellfun(@(x) strcmp(x, 'control'), evt.Modifier))
					ToggleStaticRoi(hFig);
				end
			case 'r'
				ReplicateRoi(hFig, evt);
			case 'c'
				CopyRoi(hFig);
			case 'v'
				PasteRoi(hFig);
			case 'g'
				ToggleRoiFigure(hFig);
			case 't'
				ToggleRoiDisplay(hFig);
			case 'z'
				TogglePanZoomMode(hFig);
			case 'n'
				ToggleNudgeMode(hFig);
			case 'h'
				HideExtras(hFig);
			case 'i'
				ToggleRoiInfo(hFig);
			case 'j'
				StartStopJitter(hFig);
			case 'u'
				UndoRoi(hFig);
			case 'p'
				CurrentAxes = get(hFig, 'CurrentAxes');
				selectedRoi = getappdata(hFig, 'selectedRoi');

				InterpolateRoi(CurrentAxes, selectedRoi);
% 				StartStopMovie(hFig);
			case 'a'
				InterpolateAcrossTime(hFig);
			case 'f'
% 				ToggleFreezeColormap(get(hFig, 'CurrentAxes'));
				FlipAxes(hFig);
			case 'o'
				ToggleCircle(hFig);
			case 'comma'
				KbScaleRoi(hFig, 1-prefs.roiScaleFactor);
			case 'period'
				KbScaleRoi(hFig, 1+prefs.roiScaleFactor);
			case 'k'
% 				keyboard
			case 'q'
				SetNormalMode(hFig);
				UnselectAllRois(hFig);
			case 'b'
				hColorbar = getappdata(hFig, 'hColorbar');
				if isempty(hColorbar)
					hColorbar = colorbar('XColor', [1 1 1], 'YColor', [1 1 1]);
					setappdata(hFig, 'hColorbar', hColorbar);
				else
					colorbar(hColorbar, 'off')
					setappdata(hFig, 'hColorbar', []);
				end
			case 'x'
				TogglePlotIntersect(hFig);
			case 'y'
				ToggleRoiMask(hFig);
% 				TestCase
			case 'decimal'
				MarkFrame(hFig);
			case 'f12'
				try
					prefs = getappdata(hFig, 'prefs');
					% Temporarily remove the keyboard callbacks for the editor
					kbCallbacks = prefs.kbCallbacks;
					prefs.kbCallbacks = [];
					ws = warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
					cu = onCleanup(@(x) warning(ws));

					[~, prefs] = propertiesGUI(prefs);
					prefs.kbCallbacks = kbCallbacks;
					setappdata(hFig, 'prefs', prefs);
				end
			case 'l'
				if ~iscell(prefs.transcribeRoi)
					prefs = getappdata(hFig, 'prefs');
					prefs.transcribeRoi = ~prefs.transcribeRoi;
					setappdata(hFig, 'prefs', prefs);
				end
			otherwise
%  				disp(evt.Key)
		end
	end

	% Catch key-presses in the ROI data plot window and send them to the main window
	function CatchKeyPressRoiPlot(hFig, evt)
		CatchKeyPress(getappdata(hFig, 'hMainFig'), evt);
	end

%% % % ------------------------------------------------------------------------------
function StartStopMovie(hFig)
	isPlaying = getappdata(hFig, 'isPlaying');

	if isPlaying
		setappdata(hFig, 'isPlaying', false);
		pause(1) % Poor attempt to prevent race condition
		return
	else
		setappdata(hFig, 'isPlaying', true);
		hAxis = get(hFig, 'CurrentAxes');
		for i = 1:5000  % If the stability of this loop is trusted, this could be a while(1) loop instead
			if ~ishandle(hFig)
				% Looks like someone closed the figure on us!
				return
			end
			if ~getappdata(hFig, 'isPlaying')
				break
			end
			SetTimeSlice(hFig, getappdata(hAxis, 'currentTime')+1,   getappdata(hAxis, 'currentSlice'));
			pause(1/30) % 30 fps.  Faster may be pushing it, depending on the system.
		end
	end
end

function StartStopJitter(hFig)
	n = 3; % 3 frames on either side
	isPlaying = getappdata(hFig, 'isPlaying');

	if isPlaying
		setappdata(hFig, 'isPlaying', false);
		pause(0.1) % Poor attempt to prevent race condition
		return
	else
		setappdata(hFig, 'isPlaying', true);
		hAxis = get(hFig, 'CurrentAxes');

		currentTime = getappdata(hAxis, 'currentTime');
		currentSlice = getappdata(hAxis, 'currentSlice');		
		aTimes = currentTime+[1:n n-1:-1:-n -n+1:0]; % "Yo-yo"
		for i = 1:5000 % If the stability of this loop is trusted, this could be a while(1) loop instead
			for j = 1:numel(aTimes)
				if ~ishandle(hFig)
					% Looks like someone closed the figure on us!
					return
				end
				if ~getappdata(hFig, 'isPlaying')
					SetTimeSliceImageOnly(hFig, currentTime, currentSlice);
					setappdata(hAxis, 'currentSliceImg', []);
					setappdata(hAxis, 'currentTimeImg',  []);
					break
				end
				SetTimeSliceImageOnly(hFig, aTimes(j), currentSlice);
				pause(1/30)  % 30 fps.  Faster may be pushing it, depending on the system.
			end
		end
	end
end

%% -----------------------------------------------------------------------------
function ExportImage(hFig)
	persistent lastPath

	% Dirty workaround for the fact that uigetfile() doesn't allow you to specify start location
	if exist('lastPath', 'var') && ~isempty(lastPath)
		currentDir = pwd;
		
		if ~all(lastPath == 0)
			cd(lastPath)
		end

		cu = onCleanup(@()cd(currentDir)); % Return to pwd when this function ends (even if crash/dbstop)
	end

	[fileName, pathName] = uiputfile({'*.gif','Animated GIF'; '*.jpg;*.jpeg','JPEG image'; '*.png','Portable Network Graphics file'; '*.avi', 'Motion JPEG Video'; '*.mp4', 'MPEG-4 Video'}, 'Select output file');
	lastPath = pathName;

	% Abort if dialog box was closed
	if all(fileName == 0) && all(pathName == 0)
		return
	end

	% Get relevant data
	currentSlice = getappdata(get(hFig, 'CurrentAxes'), 'currentSlice');
	currentTime  = getappdata(get(hFig, 'CurrentAxes'), 'currentTime');
	imgData      = getappdata(get(hFig, 'CurrentAxes'), 'imgData');

	[~, ~, ext] = fileparts(fileName);
	switch ext
		case '.gif'
			bHideText  = reinput('Hide text', true, 'bool');
			if numel(imgData{currentSlice}) > 1
				startFrame = reinput('Start frame', 1);
				endFrame   = reinput('End frame', numel(imgData{currentSlice}));
				skip       = reinput('Use every ''n''th frame (1 for all frames, ''inf'' for single frame)', 1);
			else
				startFrame = 1;
				endFrame   = 1;
				skip = inf;
			end

			if ~isinf(skip)
				fps        = reinput('Frame rate in fps (suggested max = 15)', 15);
				loopCount  = reinput('Number of repetitions (0 to 65536 or ''inf'')', inf);
			end
		case {'.jpg', '.jpeg'}
			bHideText  = reinput('Hide text', true, 'bool');
			if numel(imgData{currentSlice}) > 1
				startFrame = reinput('Frame number', currentTime);
			else
				startFrame = currentTime;
			end
			endFrame   = startFrame;
			skip       = inf;
		case {'.png'}
			bHideText  = reinput('Hide text', false, 'bool');
			if numel(imgData{currentSlice}) > 1
				startFrame = reinput('Frame number', currentTime);
			else
				startFrame = currentTime;
			end
			endFrame   = startFrame;
			skip       = inf;
		case {'.avi', '.mp4'}
			bHideText  = reinput('Hide text', true, 'bool');
			startFrame = reinput('Start frame', 1);
			endFrame   = reinput('End frame', numel(imgData{currentSlice}));
			skip       = reinput('Use every ''n''th frame (1 for all frames, ''inf'' for single frame)', 1);
			fps        = reinput('Frame rate in fps', 30);
	end

	if strcmpi(ext, '.avi')
		profile = 'Motion JPEG AVI'; % Compressed AVI file using Motion JPEG codec. (default)
	elseif strcmpi(ext, '.mp4')
		profile = 'MPEG-4';          % Compressed MPEG-4 file with H.264 encoding (Windows 7 and Mac OS X 10.7 only)
	end

	disp('Please wait...')

	% Hide text if currently visible
	if bHideText
		bUnhideExtras = false;
		if strcmp(get(hFig, 'Renderer'), 'zbuffer') || strcmp(get(hFig, 'Renderer'), 'opengl') % This is a more reliable indicator if everything has been hidden or not
			HideExtras(hFig);
			bUnhideExtras = true;
		end
	end

	% Explicitly hide the status text, as this is almost always unnecessary
	currStatusVis = get(getappdata(hFig, 'hStatusText'), 'Visible');
	set(getappdata(hFig, 'hStatusText'), 'Visible', 'off')

	% Explicitly hide the axes box, as this is almost always unnecessary too
	currStatusAxesBox = get(getappdata(hFig, 'hAxesBox'), 'Visible');
	set(getappdata(hFig, 'hAxesBox'), 'Visible', 'off')
	
% 	% Get relevant data
% 	currentSlice = getappdata(get(hFig, 'CurrentAxes'), 'currentSlice');
% 	currentTime  = getappdata(get(hFig, 'CurrentAxes'), 'currentTime');
% 	imgData      = getappdata(get(hFig, 'CurrentAxes'), 'imgData');
% 

	% getframe needs the figure to be on the primary screen, not on secondary monitors
	origUnitsFig = get(hFig, 'Units');
	cu1 = onCleanup(@()set(hFig, 'Units', origUnitsFig));
	
	origUnitsRoot = get(0, 'Units');
	cu2 = onCleanup(@()set(0, 'Units', origUnitsRoot));

	set(hFig, 'Units', 'pixels');
	set(0,    'Units', 'pixels');

	ScreenSize = get(0, 'ScreenSize');
	figPos     = get(hFig, 'Position');

	if figPos(1) > ScreenSize(3)
		figPos(1) = ScreenSize(3)-figPos(3);
		set(hFig, 'Position', figPos)
	end

	goodFrames = startFrame:skip:endFrame;
	frames = repmat(struct('cdata', [], 'colormap', []), [1 max(goodFrames)]);

	% capturescreen requires that the figure be fully on the main screen,
	% so move if necessary
	figPos = get(hFig, 'Position');
	figPosNew = max(figPos, [0 0 0 0]);
	if ~all(figPosNew == figPos)
		set(hFig, 'Position', figPosNew);
	end

	for iTime = goodFrames
		SetTimeSlice(hFig, iTime, currentSlice);
		drawnow; pause(0.01)
		frames(iTime) = getframe(hFig);
	end

	% Different export cases
	switch ext
		case '.gif'
			tmpData = cat(2, frames(goodFrames).cdata);

			if isequal(tmpData(:,:,1), tmpData(:,:,2)) && isequal(tmpData(:,:,1), tmpData(:,:,3))
				% Use a logical colormap for grayscale
				cMap = gray(256);
				indImg = reshape(tmpData(:,:,1), size(frames(startFrame).cdata,1), size(frames(startFrame).cdata,2), []);
			else
				% Covert from RGB to indexed
				[indImg, cMap] = rgb2ind(tmpData, 256);
				indImg = reshape(indImg, size(frames(startFrame).cdata,1), size(frames(startFrame).cdata,2), []);
			end

			if isinf(skip)
				imwrite(permute(indImg, [1 2 4 3]), cMap, fullfile(pathName, fileName), 'gif');
			else
				imwrite(permute(indImg, [1 2 4 3]), cMap, fullfile(pathName, fileName), 'gif', 'LoopCount', loopCount, 'DelayTime', 1/fps);
			end
		case {'.jpg', '.jpeg', '.png'}
			imwrite(frames(startFrame).cdata, fullfile(pathName, fileName));
		case {'.avi', '.mp4'}
			writerObj = VideoWriter(fullfile(pathName, fileName), profile);
			writerObj.FrameRate = fps;
			writerObj.Quality = 100;
			open(writerObj)
			writeVideo(writerObj,frames(goodFrames));
			close(writerObj);
	end
	disp(sprintf('Exported to "%s"', fullfile(pathName, fileName)))

	SetTimeSlice(hFig, currentTime, currentSlice);

	% Unhide text
	if bHideText && bUnhideExtras
		HideExtras(hFig);
	end

	% Restore visible status of status text
	set(getappdata(hFig, 'hStatusText'), 'Visible', currStatusVis);
	
	% Restore visible status of axes box
	set(getappdata(hFig, 'hAxesBox'), 'Visible', currStatusAxesBox);
end

%% -----------------------------------------------------------------------------
function ToggleRoiInfo(hFig)
	prefs = getappdata(hFig, 'prefs');
% 	if prefs.hideRoiLabels
% 		prefs.hideRoiLabels = false;
% 	else
% 		prefs.hideRoiLabels = true;
% 	end

	% Cycle to next mode
	switch(prefs.roiLabelMode)
		case 'none'
			prefs.roiLabelMode = 'number';
		case 'number'
			prefs.roiLabelMode = 'area';
		case 'area'
			prefs.roiLabelMode = 'signal';
		case 'signal'
			prefs.roiLabelMode = 'both';
		case 'both'
			prefs.roiLabelMode = 'none';
	end

	setappdata(hFig, 'prefs', prefs);
	UpdateAllRoiHandles(hFig)
end

function FlipAxes(hFig, evt)
	CurrentAxes = get(hFig, 'CurrentAxes');
	hAxes = getappdata(hFig, 'hAxes');
	
	currXDir = get(CurrentAxes, 'XDir');
	currYDir = get(CurrentAxes, 'YDir');

	if strcmp(currXDir, 'normal') && strcmp(currYDir, 'reverse')
		% Equivalent to fliplr
		XDir = 'reverse';
		YDir = 'reverse';
	elseif strcmp(currXDir, 'reverse') && strcmp(currYDir, 'reverse')
		% Equivalent to flipud(fliplr), or rot90(x,2)
		XDir = 'reverse';
		YDir = 'normal';
	elseif strcmp(currXDir, 'reverse') && strcmp(currYDir, 'normal')
		% Equivalent to flipud
		XDir = 'normal';
		YDir = 'normal';
	elseif strcmp(currXDir, 'normal') && strcmp(currYDir, 'normal')
		% This is the standard orientation
		XDir = 'normal';
		YDir = 'reverse';
	end

	% Set for all axes
	for hAxis = hAxes(hAxes ~= 0)'
		set(hAxis, 'XDir', XDir)
		set(hAxis, 'YDir', YDir)
	end
end

function InterpolateImages(hFig)
	CurrentAxes = get(hFig, 'CurrentAxes');
	hAxes       = getappdata(hFig,        'hAxes');
	currentTime = getappdata(CurrentAxes, 'currentTime');
	
	for hAxis = reshape(hAxes(hAxes ~= 0), 1, [])
		imgData      = getappdata(hAxis, 'imgData');
		for iSli = 1:numel(imgData)
			for iTime = 1:numel(imgData{iSli})
				if ~isempty(imgData{iSli}{iTime})
					imgData{iSli}{iTime} = imresize(imgData{iSli}{iTime}, 2);
				end
			end
		end
		setappdata(hAxis, 'imgData', imgData);
		setappdata(hAxis, 'currentTime', 0); % Dummy value to force a reload of the image
	end

	SetTimeSlice(hFig, currentTime, getappdata(CurrentAxes, 'currentSlice'));
end

%% -----------------------------------------------------------------------------
function MouseScroll(hFig, evt)
	% In R2014b+ evt.VerticalScrollCount is read-only
	VerticalScrollCount = evt.VerticalScrollCount;

	% Apple Magic Mouse and Magic Trackpad seem to always return 2 scrolls at a
	% time. Programmatically limit this to 1.
	if VerticalScrollCount < 0
		VerticalScrollCount = -1;
	else
		VerticalScrollCount = 1;
	end
	
	hAxis = get(hFig, 'CurrentAxes');
	CurrentModifier = get(hFig, 'CurrentModifier');
	
	if any(cellfun(@(x) strcmp(x, 'shift'), CurrentModifier))
		CurrentAxes = get(hFig, 'CurrentAxes');
		szImg = size(get(findobj(get(CurrentAxes, 'Children'), 'Type', 'image'), 'CData'));
		szAxis = [diff(get(CurrentAxes, 'YLim')) diff(get(CurrentAxes, 'XLim'))];
		currentZoom = max(szImg ./ szAxis);

		setappdata(hFig, 'InitPoint', get(CurrentAxes, 'CurrentPoint'))
		setappdata(hFig, 'InitZoom',  currentZoom)

		AdjustZoom(hFig, [], 1 - VerticalScrollCount/10);

	else
		SetTimeSlice(hFig, getappdata(hAxis, 'currentTime')-VerticalScrollCount, getappdata(hAxis, 'currentSlice'));
	end
end


%% ------------------------------------------------------------
function MouseDownSelect(hFig, evt)
% Called on any mouse-down event
% - Update selected subplot if applicable
% - Update selected ROI if applicable

	CurrentAxes     = get(hFig,        'CurrentAxes');
	CurrentObject   = get(hFig,        'CurrentObject');
	lastCurrentAxes = getappdata(hFig, 'lastCurrentAxes');
	
	SelectSubplot(hFig, CurrentAxes);
	SelectRoi(hFig, lastCurrentAxes, CurrentAxes, CurrentObject);
end




% --- 1. Show which axes is selected with a red box --------------------------
function SelectSubplot(hFig, CurrentAxes)
	lastCurrentAxes = getappdata(hFig, 'lastCurrentAxes');
	hAxesBox        = getappdata(hFig, 'hAxesBox');

	if (CurrentAxes == lastCurrentAxes) || isempty(hAxesBox)
		return
	end

	setappdata(hFig, 'lastCurrentAxes', CurrentAxes);

	axesPosition = get(CurrentAxes, 'Position');
	set(hAxesBox, 'XData', axesPosition(1) + [0 axesPosition(3) axesPosition(3) 0 0] + 0.001*[1 -1 -1 1 1])
	set(hAxesBox, 'YData', axesPosition(2) + [0 0 axesPosition(4) axesPosition(4) 0] + 0.001*[1 1 -1 -1 1])

	% De-select all ROIs
	SelectRoi(hFig, lastCurrentAxes, CurrentAxes, []);
end

% --- 2. Highlights/unhighlights appropriate ROIs ----------------------------
function SelectRoi(hFig, lastCurrentAxes, CurrentAxes, CurrentObject)
	roiHandles     = getappdata(CurrentAxes, 'roiHandles');
	oldSelectedRoi = getappdata(hFig,        'selectedRoi');
	newSelectedRoi = [];

	% Find (and highlight) newly selected ROI
	for iRoi = 1:length(roiHandles)
		tmpHandles = [roiHandles(iRoi).line roiHandles(iRoi).vertices roiHandles(iRoi).center roiHandles(iRoi).label];

		if isempty(tmpHandles)
			continue
		end
		
		% This no longer works in R2014b because handles aren't numeric. New code
		% should be ok as long as CurrentObject isn't an array
% 		if intersect(tmpHandles, CurrentObject)
		if any(arrayfun(@(x) any(tmpHandles == x), CurrentObject))
			set(roiHandles(iRoi).line,     'LineStyle', '--');
			set(roiHandles(iRoi).vertices, 'Marker',   '+');
			
			if ismac
				set(roiHandles(iRoi).vertices, 'MarkerSize', 6)
			end
			newSelectedRoi = [newSelectedRoi iRoi];
		end
	end

	% Unhighlight previously selected ROI
	if (lastCurrentAxes == CurrentAxes) && ...
	   (numel(oldSelectedRoi) == numel(newSelectedRoi)) && ...
	   (numel(intersect(oldSelectedRoi, newSelectedRoi)) == numel(newSelectedRoi))
		% Same ROI on the same subplot; do nothing
	elseif isempty(oldSelectedRoi)
		% No ROI to unhighlight
	else
		unhighlightRois = setdiff(oldSelectedRoi, newSelectedRoi);
		oldRoiHandles = getappdata(lastCurrentAxes, 'roiHandles');
		unhighlightRois(unhighlightRois > numel(oldRoiHandles)) = [];
		if ~isempty(unhighlightRois)
% 			oldRoiHandles = getappdata(lastCurrentAxes, 'roiHandles');
% 			unhighlightRois(unhighlightRois > numel(oldRoiHandles)) = [];
			set([oldRoiHandles(unhighlightRois).line],     'LineStyle', '-');
			set([oldRoiHandles(unhighlightRois).vertices], 'Marker',    '.');

			if ismac
				set([oldRoiHandles(unhighlightRois).vertices], 'MarkerSize', 12)
			end
		end
	end

	setappdata(hFig, 'selectedRoi', newSelectedRoi);
end

function UnselectAllRois(hFig)
	hAxes = getappdata(hFig, 'hAxes');
% 	CurrentAxes = get(hFig, 'CurrentAxes');

	for jAxis = reshape(hAxes(hAxes ~= 0), 1, [])
		roiHandles  = getappdata(jAxis, 'roiHandles');

		for iRoi = 1:numel(roiHandles)
			set([roiHandles(iRoi).line],     'LineStyle', '-');
			set([roiHandles(iRoi).vertices], 'Marker',    '.');

			if ismac
				set([roiHandles(iRoi).vertices], 'MarkerSize', 12)
			end
		end
	end

	setappdata(hFig, 'selectedRoi', []);
end

	% ------------------------------------------------------------
	function DeleteRois(hFig, aDeleteRois)
	% - Deletes ROI data for 'aRois'
	% - Sets the appropriate ROI handles to display nothing
	% - Sets 'selectedRoi' to empty

		CurrentAxes = get(hFig, 'CurrentAxes');
		staticRoi   = getappdata(hFig, 'staticRoi');
		prefs       = getappdata(hFig, 'prefs');
		hAxes       = getappdata(hFig, 'hAxes');

		% Reset selectedRoi as needed
		selectedRoi = getappdata(hFig, 'selectedRoi');
		if isempty(aDeleteRois) || any(aDeleteRois == selectedRoi)
			setappdata(hFig, 'selectedRoi', []);
		end

		% transcribeRoi dictates which axes we act on
		if iscell(prefs.transcribeRoi)
			currentSubplot = find(hAxes == CurrentAxes);
			goodInds = cellfun(@(x) any(x == currentSubplot), prefs.transcribeRoi);
			if any(goodInds)
				aSubplots = prefs.transcribeRoi{goodInds};
				aAxes     = hAxes(aSubplots);
			else
				aAxes = CurrentAxes;
			end
		elseif prefs.transcribeRoi
			aAxes = hAxes(hAxes ~= 0);
		elseif ~prefs.transcribeRoi
			aAxes = CurrentAxes;
		end

		aAxes = reshape(aAxes, 1, []);
		
		for jAxis = aAxes
			roi          = getappdata(jAxis, 'roi');
			roiMask      = getappdata(jAxis, 'roiMask');
			roiData      = getappdata(jAxis, 'roiData');
			roiHandles   = getappdata(jAxis, 'roiHandles');
			currentSlice = getappdata(jAxis, 'currentSlice');
			currentTime  = getappdata(jAxis, 'currentTime');

			% Save the current ROI setfor undo
			setappdata(jAxis, 'prevRoi', roi);
			
			if isempty(aDeleteRois)
				aRois = find(cellfun(@(x) ~isempty(x), roi{currentSlice}{currentTime}));
			else
				aRois = aDeleteRois;
			end

			% Change the ROI handles to display nothing
			for iRoi = aRois
				set(roiHandles(iRoi).line,     'XData', [],   'YData', []);
				set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
				set(roiHandles(iRoi).center,   'XData', [],   'YData', []);
				if ~strcmp(prefs.roiLabelMode, 'none')
					set(roiHandles(iRoi).label,    'Position', [0 0 0]);
				end

				% Unhighlight the ROI
				set(roiHandles(iRoi).line,     'LineStyle', '-');
				set(roiHandles(iRoi).vertices, 'Marker',    '.');

				% Clear the appropriate data
				roi{    currentSlice}{currentTime}{iRoi} = [];
				roiMask{currentSlice}{currentTime}{iRoi} = [];
				roiData{currentSlice}{currentTime}{iRoi} = [];

				% For staticRoi, delete it on all other time frames
				if staticRoi
					for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
						roi{    currentSlice}{iTime}{iRoi} = [];
						roiMask{currentSlice}{iTime}{iRoi} = [];
						roiData{currentSlice}{iTime}{iRoi} = [];
					end
				end
			end

			setappdata(jAxis, 'roi',        roi);
			setappdata(jAxis, 'roiMask',    roiMask);
			setappdata(jAxis, 'roiData',    roiData);
		end

		% Update the ROI data plot
		UpdateRoiFigure(hFig);

		% Apply to SimpleViewer3D, if applicable
		hSV3D = getappdata(hFig, 'hSV3D');
		if ~isempty(hSV3D) && ishandle(hSV3D)
			roi = getappdata(hAxes(1), 'roi');
			setappdata(hSV3D, 'roi', roi);
			fcnHandles = getappdata(hSV3D, 'fcnHandles');
			for jRoi = aDeleteRois
				fcnHandles.UpdateRoiHandles(hSV3D, jRoi);
			end
		end
	end

% 	function DeleteAllRois(hFig)
% 	% - Deletes ROI data for selected ROI
% 	% - Sets the appropriate ROI handles to display nothing
% 	% - Sets 'selectedRoi' to empty
% 
% 		CurrentAxes = get(hFig, 'CurrentAxes');
% 		staticRoi   = getappdata(hFig, 'staticRoi');
% 		prefs       = getappdata(hFig, 'prefs');
% 		hAxes       = getappdata(hFig, 'hAxes');
% 
% 		for jAxis = hAxes(:)'
% 			if (jAxis == 0)
% 				continue
% 			end
% 
% 			% Don't do other subplots if transcribeRoi is off
% 			if ~prefs.transcribeRoi && (jAxis ~= CurrentAxes)
% 				continue
% 			end
% 
% 			roi          = getappdata(jAxis, 'roi');
% 			roiMask      = getappdata(jAxis, 'roiMask');
% 			roiData      = getappdata(jAxis, 'roiData');
% 			roiHandles   = getappdata(jAxis, 'roiHandles');
% 			currentSlice = getappdata(jAxis, 'currentSlice');
% 			currentTime  = getappdata(jAxis, 'currentTime');
% 
% 			% Change the ROI handles to display nothing
% 			for iRoi = find(cellfun(@(x) ~isempty(x), roi{currentSlice}{currentTime}))
% 				set(roiHandles(iRoi).line,     'XData', [],   'YData', []);
% 				set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
% 				set(roiHandles(iRoi).center,   'XData', [],   'YData', []);
% 				if ~prefs.hideRoiLabels
% 					set(roiHandles(iRoi).label,    'Position', [0 0 0]);
% 				end
% 
% 				% Unhighlight the ROI
% 				set(roiHandles(iRoi).line,     'LineStyle', '-');
% 				set(roiHandles(iRoi).vertices, 'Marker',    '.');
% 
% 				% Clear the appropriate data
% 				roi{    currentSlice}{currentTime}{iRoi} = [];
% 				roiMask{currentSlice}{currentTime}{iRoi} = [];
% 				roiData{currentSlice}{currentTime}{iRoi} = [];
% 
% 				% For staticRoi, delete it on all other time frames
% 				if staticRoi
% 					for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 						roi{    currentSlice}{iTime}{iRoi} = [];
% 						roiMask{currentSlice}{iTime}{iRoi} = [];
% 						roiData{currentSlice}{iTime}{iRoi} = [];
% 					end
% 				end
% 			end
% 
% 			setappdata(jAxis, 'roi',        roi);
% 			setappdata(jAxis, 'roiMask',    roiMask);
% 			setappdata(jAxis, 'roiData',    roiData);
% 		end
% 		setappdata(hFig, 'selectedRoi', []);
% 
% 		% Update the ROI data plot
% 		UpdateRoiFigure(hFig);
% 	end
% 
% 	% ------------------------------------------------------------
% 	function DeleteSelectedRoi(hFig)
% 	% - Deletes ROI data for selected ROI
% 	% - Sets the appropriate ROI handles to display nothing
% 	% - Sets 'selectedRoi' to empty
% 
% 		hAxis       = get(hFig, 'CurrentAxes');
% 		selectedRoi = getappdata(hFig, 'selectedRoi');
% 		staticRoi   = getappdata(hFig, 'staticRoi');
% 		prefs       = getappdata(hFig, 'prefs');
% 
% 		if isempty(selectedRoi) || ~isempty(getappdata(hAxis, 'inProgressRoi'))
% 			return
% 		end
% 
% 		roi          = getappdata(hAxis, 'roi');
% 		roiMask      = getappdata(hAxis, 'roiMask');
% 		roiData      = getappdata(hAxis, 'roiData');
% 		roiHandles   = getappdata(hAxis, 'roiHandles');
% 		currentSlice = getappdata(hAxis, 'currentSlice');
% 		currentTime  = getappdata(hAxis, 'currentTime');
% 
% 		for jRoi = selectedRoi
% 			% Change the ROI handles to display nothing
% 			set(roiHandles(jRoi).line,     'XData', [],   'YData', []);
% 			set(roiHandles(jRoi).vertices, 'XData', [],   'YData', []);
% 			set(roiHandles(jRoi).center,   'XData', [],   'YData', []);
% 			if ~prefs.hideRoiLabels
% 				set(roiHandles(jRoi).label,  'Position', [0 0 0]);
% 			end
% 
% 			% Unhighlight the ROI
% 			set(roiHandles(jRoi).line,     'LineStyle', '-');
% 			set(roiHandles(jRoi).vertices, 'Marker',   '.');
% 
% 			% Clear the appropriate data
% 			roi{    currentSlice}{currentTime}{jRoi} = [];
% 			roiMask{currentSlice}{currentTime}{jRoi} = [];
% 			roiData{currentSlice}{currentTime}{jRoi} = [];
% 
% 			% For staticRoi, delete it on all other time frames
% 			if staticRoi
% 				for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 					roi{    currentSlice}{iTime}{jRoi} = [];
% 					roiMask{currentSlice}{iTime}{jRoi} = [];
% 					roiData{currentSlice}{iTime}{jRoi} = [];
% 				end
% 			end
% 		end
% 
% 		setappdata(hAxis, 'roi',        roi);
% 		setappdata(hAxis, 'roiMask',    roiMask);
% 		setappdata(hAxis, 'roiData',    roiData);
% 		setappdata(hFig, 'selectedRoi', []);
% 
% 		% Update the ROI data plot
% 		UpdateRoiFigure(hFig);
% 
% % % 		% Apply to another SimpleViewer, if applicable
% % % 		hSV2D = getappdata(hFig, 'hSV2D');
% % % 		if ~isempty(hSV2D) && ishandle(hSV2D)
% % % 			
% % % % 			hAxes = getappdata(hFig, 'hAxes');
% % % % 			if (CurrentAxes == hAxes(1))
% % % 				setappdata(hSV2D, 'roi', roi);
% % % 				fcnHandles = getappdata(hSV2D, 'fcnHandles');
% % % 				for jRoi = selectedRoi
% % % 					fcnHandles.UpdateRoiHandles(hSV2D, jRoi);
% % % 				end
% % % % 			end
% % % 		end
% 
% 		% Apply to SimpleViewer3D, if applicable
% 		hSV3D = getappdata(hFig, 'hSV3D');
% 		if ~isempty(hSV3D) && ishandle(hSV3D)
% % 			hAxes = getappdata(hFig, 'hAxes');
% % 			if (CurrentAxes == hAxes(1))
% 				setappdata(hSV3D, 'roi', roi);
% 				fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 				for jRoi = selectedRoi
% 					fcnHandles.UpdateRoiHandles(hSV3D, jRoi);
% 				end
% % 			end
% 		end
% 	end

%% -----------------------------------------------------------------------------
function DoStaticAndTranscribeAnd3D(hFig, selectedRois, bDontUpdateData)
% 	% Propagate (if static) and copy (if transcribe) ROIs
	% - Requires roi be properly set for CurrentAxes and currentSlice/currentTime

	CurrentAxes     = get(hFig, 'CurrentAxes');
	staticRoi       = getappdata(hFig, 'staticRoi');
	prefs           = getappdata(hFig, 'prefs');
	hAxes           = getappdata(hFig, 'hAxes');

	currentSlice    = getappdata(CurrentAxes, 'currentSlice');
	currentTime     = getappdata(CurrentAxes, 'currentTime');
	roi             = getappdata(CurrentAxes, 'roi');

	if ~exist('bDontUpdateData', 'var')
		bDontUpdateData = false;
	end

	% transcribeRoi can be:
	% - false:       No transcribing
	% - true:        Transcribe to all plots
	% - cell array:  Each cell is an array containing subplot indices.  Link each subplot in each array

	if iscell(prefs.transcribeRoi)
		currentSubplot = find(hAxes == CurrentAxes);
		goodInds = cellfun(@(x) any(x == currentSubplot), prefs.transcribeRoi);
		if any(goodInds)
			aSubplots = prefs.transcribeRoi{goodInds};
			aAxes     = hAxes(aSubplots);
		else
			aAxes = CurrentAxes;
		end
	elseif prefs.transcribeRoi
		aAxes = hAxes(hAxes ~= 0);
	elseif ~prefs.transcribeRoi
		aAxes = CurrentAxes;
	end
	
	aAxes = reshape(aAxes, 1, []);

% 	for jAxis = hAxes(:)'
% 		if (jAxis == 0)
% 			continue
% 		end
% 
% 		% Don't do other subplots if transcribeRoi is off
% 		if ~prefs.transcribeRoi && (jAxis ~= CurrentAxes)
% 			continue
% 		end

	for jAxis = aAxes
		tRoi = getappdata(jAxis, 'roi');

		% Range check (these are clipped, not cyclic)
		tSlice = min([currentSlice length(tRoi        )]);
		tTime  = min([currentTime  length(tRoi{tSlice})]);

		tRoi{tSlice}{tTime}(selectedRois) = roi{currentSlice}{currentTime}(selectedRois);

		% For staticRoi, propagate it to all other time frames
		if staticRoi
			for iTime = [1:tTime-1 tTime+1:numel(tRoi{tSlice})]
				tRoi{tSlice}{iTime}(selectedRois) = tRoi{tSlice}{tTime}(selectedRois);
			end
		end

		setappdata(jAxis, 'roi', tRoi);
		for iRoi = selectedRois
			if bDontUpdateData
				UpdateRoiHandles(jAxis, iRoi);
			else
				if staticRoi
					UpdateRoiData(jAxis, iRoi);
				else
					UpdateRoiData(jAxis, iRoi, tTime);
				end
			end
		end
	end

	% Apply to SimpleViewer3D, if applicable
	hSV3D = getappdata(hFig, 'hSV3D');
	if ~isempty(hSV3D) && ishandle(hSV3D)
		if (CurrentAxes == hAxes(1))
			setappdata(hSV3D, 'roi', getappdata(CurrentAxes, 'roi'));
			fcnHandles = getappdata(hSV3D, 'fcnHandles');
			% FIXME: These functions should really support multiple selectedRois...
			for iRoi = selectedRois
				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
			end
		end
	end

% % % 	% For staticRoi, propagate it to all other time frames
% % % 	if staticRoi
% % % 		for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% % % 			roi{currentSlice}{iTime}(selectedRoi) = roi{currentSlice}{currentTime}(selectedRoi);
% % % 		end
% % % 	end
% % % 	setappdata(CurrentAxes, 'roi', roi);
% % % 
% % % 	% Transcribe to other subplots
% % % 	if prefs.transcribeRoi
% % % 		hAxes = getappdata(hFig, 'hAxes');
% % % 		for jAxis = hAxes(:)'
% % % 			if (jAxis == 0) || (jAxis == CurrentAxes)
% % % 				continue
% % % 			end
% % % 			tRoi = getappdata(jAxis, 'roi');
% % % 			tRoi{currentSlice}{currentTime}{selectedRoi} = roi{currentSlice}{currentTime}{selectedRoi};
% % % 
% % % 			% For staticRoi, propagate it to all other time frames
% % % 			if staticRoi
% % % 				for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% % % 					tRoi{currentSlice}{iTime}{selectedRoi} = tRoi{currentSlice}{currentTime}{selectedRoi};
% % % 				end
% % % 			end
% % % 
% % % 			setappdata(jAxis, 'roi', tRoi);
% % % 			UpdateRoiHandles(jAxis, selectedRoi);
% % % 		end
% % % 	end
end

% ------------------------------------------------------------------------------
function KbMoveRoi(hFig, xTrans, yTrans)
% Move a selected ROI via keyboard

	CurrentAxes     = get(hFig, 'CurrentAxes');
	selectedRoi     = getappdata(hFig, 'selectedRoi');
	currentSlice    = getappdata(CurrentAxes, 'currentSlice');
	currentTime     = getappdata(CurrentAxes, 'currentTime');
	roi             = getappdata(CurrentAxes, 'roi');

	if isempty(selectedRoi)
		return
	end

	for iRoi = selectedRoi
		roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi} + repmat([xTrans yTrans], [size(roi{currentSlice}{currentTime}{iRoi},1) 1]);
	end
	setappdata(CurrentAxes, 'roi', roi);

	% Propagate and transcribe moved ROI as necessary
	DoStaticAndTranscribeAnd3D(hFig, selectedRoi);

% 	% TODO: Do we still need to do this?
% 	% FIXME: These functions should really support multiple selectedRois...
% 	for iRoi = selectedRoi
% 		UpdateRoiData(hFig, iRoi);
% 	end
% 
% 	% Apply to SimpleViewer3D, if applicable
% 	hSV3D = getappdata(hFig, 'hSV3D');
% 	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		hAxes = getappdata(hFig, 'hAxes');
% 		if (CurrentAxes == hAxes(1))
% 			setappdata(hSV3D, 'roi', roi);
% 			fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 			% FIXME: These functions should really support multiple selectedRois...
% 			for iRoi = selectedRoi
% 				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 			end
% 		end
% 	end
	return
end

% ------------------------------------------------------------------------------
function KbScaleRoi(hFig, scale)
% Move a selected ROI via keyboard

	CurrentAxes     = get(hFig, 'CurrentAxes');
	selectedRoi     = getappdata(hFig, 'selectedRoi');
	currentSlice    = getappdata(CurrentAxes, 'currentSlice');
	currentTime     = getappdata(CurrentAxes, 'currentTime');
	roi             = getappdata(CurrentAxes, 'roi');

	if isempty(selectedRoi)
		return
	end

	% Scale ROI roi relative to the center
	for iRoi = selectedRoi
		roiVertices = roi{currentSlice}{currentTime}{iRoi};
		center = mean(roiVertices,1);
		roiVertices = (roiVertices - repmat(center, [size(roiVertices,1) 1])) * scale + repmat(center, [size(roiVertices,1) 1]);
		roi{currentSlice}{currentTime}{iRoi} = roiVertices;
	end
	setappdata(CurrentAxes, 'roi', roi);

	% Propagate and transcribe moved ROI as necessary
	DoStaticAndTranscribeAnd3D(hFig, selectedRoi);

% 	% FIXME: These functions should really support multiple selectedRois...
% 	for iRoi = selectedRoi
% 		UpdateRoiData(hFig, iRoi);
% 	end
% 
% 	% Apply to SimpleViewer3D, if applicable
% 	hSV3D = getappdata(hFig, 'hSV3D');
% 	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		hAxes = getappdata(hFig, 'hAxes');
% 		if (CurrentAxes == hAxes(1))
% 			setappdata(hSV3D, 'roi', roi);
% 			fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 			% FIXME: These functions should really support multiple selectedRois...
% 			for iRoi = selectedRoi
% 				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 			end
% 		end
% 	end

	return
end

	%% ------------------------------------------------------------
	function CopyRoi(hFig, useGlobal)
	% - Copy selected ROI into clipboard (if appropriate)

		CurrentAxes  = get(hFig, 'CurrentAxes');
		selectedRoi  = getappdata(hFig, 'selectedRoi');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');

		if isempty(selectedRoi)
			return
		end

		roi = getappdata(CurrentAxes, 'roi');
		
		% This is necessary because we allow an ROI number to be selected, even if
		% it hasn't yet been drawn
		if selectedRoi > numel(roi{currentSlice}{currentTime})
			return
		end

		if exist('useGlobal', 'var') && useGlobal
			global gClipboardRoi  %#ok<TLEV>  GLOBAL could be very inefficient unless it is a top-level statement in its function.
			gClipboardRoi = roi{currentSlice}{currentTime}(selectedRoi);
		else
			setappdata(hFig, 'clipboardRoi', roi{currentSlice}{currentTime}(selectedRoi));
		end
	end

	% ------------------------------------------------------------
	function PasteRoi(hFig, useGlobal)
	% - Paste ROI from clipboard (if appropriate)

		CurrentAxes  = get(hFig, 'CurrentAxes');
		roi          = getappdata(CurrentAxes, 'roi');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');
		staticRoi    = getappdata(hFig,  'staticRoi');

		if exist('useGlobal', 'var') && useGlobal
			global gClipboardRoi  %#ok<TLEV>  GLOBAL could be very inefficient unless it is a top-level statement in its function.
			clipboardRoi = gClipboardRoi;
		else
			clipboardRoi = getappdata(hFig, 'clipboardRoi');
		end

		% Abort if the clipboard is empty
		if isempty(clipboardRoi)
			return
		end
		
		aRoi = zeros(1, numel(clipboardRoi));
		for iRoi = 1:numel(clipboardRoi)
			% Next ROI is the first empty slot (i.e. if ROI 2 of 3 is deleted, the next ROI drawn is index 2)
			jRoi = find(cellfun(@(x) isempty(x), roi{currentSlice}{currentTime}), 1);
			if isempty(jRoi)
				jRoi = length(roi{currentSlice}{currentTime})+1;
			end

			roi{currentSlice}{currentTime}(jRoi) = clipboardRoi(iRoi);
			
			% Keep track of the ROIs we've changed so we can update its data and propagate if necessary
			aRoi(iRoi) = jRoi;
		end
		setappdata(CurrentAxes, 'roi', roi);
		DoStaticAndTranscribeAnd3D(hFig, aRoi)

		% Select the ROI too
		roiHandles     = getappdata(CurrentAxes, 'roiHandles');
		SelectRoi(hFig, CurrentAxes, CurrentAxes, roiHandles(aRoi(end)).center)

% 		aRoi = zeros(1, numel(clipboardRoi));
% 		for iRoi = 1:numel(clipboardRoi)
% 			% Next ROI is the first empty slot (i.e. if ROI 2 of 3 is deleted, the next ROI drawn is index 2)
% 			jRoi = find(cellfun(@(x) isempty(x), roi{currentSlice}{currentTime}), 1);
% 			if isempty(jRoi)
% 				jRoi = length(roi{currentSlice}{currentTime})+1;
% 			end
% 
% 			% Placeholders to allow loop to work
% 			roi{currentSlice}{currentTime}{jRoi} = 1;
% 			aRoi(iRoi) = jRoi;
% 		end
% 
% 		% The only thing we need to add is 'roi'. UpdateRoiData will take care of the rest
% 		if staticRoi
% 			for iTime = 1:length(roi{currentSlice})
% 				roi{currentSlice}{iTime}(aRoi) = clipboardRoi;
% 			end
% 		else
% 			roi{currentSlice}{currentTime}(aRoi) = clipboardRoi;
% 		end
% 
% 		setappdata(CurrentAxes, 'roi', roi);
% 		for jRoi = aRoi
% 			UpdateRoiData(hFig, jRoi);
% 		end
	end

	function ReplicateRoi(hFig, evt)
	% - Replicates selected ROI to all other time frames (all ROIs if none are selected)
		CurrentAxes  = get(hFig, 'CurrentAxes');
		selectedRoi  = getappdata(hFig, 'selectedRoi');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');
		roi          = getappdata(CurrentAxes, 'roi');

		% Get the ROI(s) to copy
		if ~isempty(selectedRoi)
			targetInds = selectedRoi;
		else
			targetInds = find(cellfun(@(x) ~isempty(x), roi{currentSlice}{currentTime}));
			if isempty(targetInds)
				return
			end
		end
		targetRois = roi{currentSlice}{currentTime}(targetInds);

		% Actually replicate now
		for iTime = [1:currentTime-1 currentTime+1:numel(roi{currentSlice})]
			roi{currentSlice}{iTime}(targetInds) = targetRois;
		end

		setappdata(CurrentAxes, 'roi', roi);
		
		% This call is wasteful in that it staticRoi re-does the replication.
		% However, this is a clean way of handling transcribeRoi and hSV3D
		DoStaticAndTranscribeAnd3D(hFig, targetInds)

% 		for jRoi = targetInds
% 			UpdateRoiData(hFig, jRoi);
% 		end
	end

	
	function ToggleRoiMask(hFig)
		% Use ROIs to mask image data.  If <5 ROIs, assume they're drawn for
		% short-axis myocardial analysis

		hAxes = getappdata(hFig, 'hAxes');
		prefs = getappdata(hFig, 'prefs');

		for jAxis = reshape(hAxes(hAxes ~= 0), 1, [])
			currentSlice     = getappdata(jAxis, 'currentSlice');
			currentTime      = getappdata(jAxis, 'currentTime');
			imgData          = getappdata(jAxis, 'imgData');
			nonMaskedImgData = getappdata(jAxis, 'nonMaskedImgData');
			roi              = getappdata(jAxis, 'roi');
			tmpRois          = roi{currentSlice}{currentTime};

			if ~isempty(nonMaskedImgData)
				% Already showing masked image, so put the stored images back
				setappdata(jAxis, 'imgData',          nonMaskedImgData);
				setappdata(jAxis, 'nonMaskedImgData', []);
			else
				% Create mask based on ROIs
				nonMaskedImgData = imgData;
				for iTime = 1:numel(imgData{currentSlice})
					if isempty(imgData{currentSlice}{iTime}) || isempty(roi{currentSlice}{iTime})
						continue
					end

					tmpImg  = imgData{currentSlice}{iTime};
					tmpRois = roi{currentSlice}{iTime};
					sz      = size(tmpImg);

					% Re-cast tmpImg to support NaNs if necessary
					if ~any(strcmp({'double', 'single'}, class(tmpImg)))
						tmpImg = single(tmpImg);
					end

					% Remove blank ROIs
					indsBad = cellfun(@(x) isempty(x), tmpRois);
					tmpRois(indsBad) = [];

					% Remove line segments
					indsBad = cellfun(@(x) all(isnan(x(end,:))), tmpRois);
					tmpRois(indsBad) = [];

					% Remove points and lines
					indsBad = cellfun(@(x) size(x,1) <= 2, tmpRois);
					tmpRois(indsBad) = [];

					% If all the areas are the same, then it's probably not SAX ROIs
					areas = cellfun(@(x) polyarea(x(:,1), x(:,2)), tmpRois);
					if all(abs(diff(areas)) < 1e-6)
						prefs.assumeMyoSaxROIs = false;
					end

					if prefs.assumeMyoSaxROIs && exist('IdentifySaxRois.m', 'file')
						if (numel(tmpRois) == 2)
							tmpRois{end+1} = nan(1,2);  % Add a dummy RVI point for IdentifySaxRois
						end

						if (numel(tmpRois) == 3) || (numel(tmpRois) == 4)
							% 3 ROIs is epi, endo, RVI (or epi, endo, blood)
							% 4 ROIs is epi, endo, RVI, blood pool
							inds = IdentifySaxRois(tmpRois);

							iEpi  = inds(1);
							iEndo = inds(2);

							mask = (poly2mask(tmpRois{iEpi }(:,1), tmpRois{iEpi }(:,2), sz(1), sz(2)) - ...
											poly2mask(tmpRois{iEndo}(:,1), tmpRois{iEndo}(:,2), sz(1), sz(2))) > 0;

							if numel(inds) > 3
								% Blood ROI present
								iBlood = inds(4);
								mask = (mask + poly2mask(tmpRois{iBlood}(:,1), tmpRois{iBlood}(:,2), sz(1), sz(2))) > 0;
							end
						else
							% Arbitrary number of ROIs -- add them all up
							cMask = cellfun(@(x) poly2mask(x(:,1), x(:,2), sz(1), sz(2)), tmpRois, 'UniformOutput', false);
							mask = (sum(cat(3, cMask{:}),3)) > 0;
						end
					else
						% Arbitrary number of ROIs -- add them all up
						cMask = cellfun(@(x) poly2mask(x(:,1), x(:,2), sz(1), sz(2)), tmpRois, 'UniformOutput', false);
						mask = (sum(cat(3, cMask{:}),3)) > 0;
					end

					tmpImg(~mask) = nan;
					imgData{currentSlice}{iTime} = tmpImg;
				end

				setappdata(jAxis, 'nonMaskedImgData', nonMaskedImgData);
				setappdata(jAxis, 'imgData',          imgData);
			end

			for iRoi = 1:numel(tmpRois)
				UpdateRoiData(jAxis, iRoi);
			end
			setappdata(jAxis, 'currentTime',    0); % Dummy value to force a reload of the image
		end
		% Need to refresh images
		SetTimeSlice(hFig, currentTime, currentSlice);
	end

% ------------------------------------------------------------------------------
% This version only masks the current time
% 	function ToggleRoiMask(hFig)
% 		% Use ROIs to mask image data.  If <5 ROIs, assume they're drawn for
% 		% short-axis myocardial analysis
% 
% 		hAxes = getappdata(hFig, 'hAxes');
% 
% 		for jAxis = reshape(hAxes(hAxes ~= 0), 1, [])
% 			currentSlice   = getappdata(jAxis, 'currentSlice');
% 			currentTime    = getappdata(jAxis, 'currentTime');
% 			imgData        = getappdata(jAxis, 'imgData');
% 			nonMaskedImage = getappdata(jAxis, 'nonMaskedImage');
% 			roi            = getappdata(jAxis, 'roi');
% 			tmpRois        = roi{currentSlice}{currentTime};
% 
% 			if ~isempty(nonMaskedImage)
% 				% Already showing masked image, so put the stored images back
% 				imgData{currentSlice}{currentTime} = nonMaskedImage;
% 				setappdata(jAxis, 'nonMaskedImage', []);
% 			else
% 				% Create mask based on ROIs
% 				nonMaskedImage = imgData{currentSlice}{currentTime};
% 				sz             = size(nonMaskedImage);
% 				maskedImage    = nonMaskedImage;
% 
% 				% Identify the ROIs
% 				switch numel(tmpRois)
% 					case 1
% 						% 1 ROI, so just mask that
% 						mask = poly2mask(tmpRois{1}(:,1), tmpRois{1}(:,2), sz(1), sz(2));
% 
% 					case 2
% 						% 2 ROIs is endo and epicardium
% 						areas = cellfun(@(x) polyarea(x(:,1), x(:,2)), tmpRois);
% 						iEpi  = find(areas == max(areas), 1);
% 						iEndo = find(areas == min(areas), 1);
% 
% 						mask = (poly2mask(tmpRois{iEpi }(:,1), tmpRois{iEpi }(:,2), sz(1), sz(2)) - ...
% 						        poly2mask(tmpRois{iEndo}(:,1), tmpRois{iEndo}(:,2), sz(1), sz(2))) > 0;
% 
% 					case 3
% 						% 3 ROIs is epi, endo, RVI
% 						inds = IdentifySaxRois(tmpRois);
% 						iEpi  = inds(1);
% 						iEndo = inds(2);
% 
% 						mask = (poly2mask(tmpRois{iEpi }(:,1), tmpRois{iEpi }(:,2), sz(1), sz(2)) - ...
% 						        poly2mask(tmpRois{iEndo}(:,1), tmpRois{iEndo}(:,2), sz(1), sz(2))) > 0;
% 
% 					case 4
% 						% 4 ROIs is epi, endo, RVI, blood pool
% 						inds = IdentifySaxRois(tmpRois);
% 						iEpi   = inds(1);
% 						iEndo  = inds(2);
% 						iBlood = inds(4);
% 
% 						mask = (poly2mask(tmpRois{iEpi }(:,1), tmpRois{iEpi }(:,2), sz(1), sz(2)) - ...
% 						        poly2mask(tmpRois{iEndo}(:,1), tmpRois{iEndo}(:,2), sz(1), sz(2))) > 0;
% 
% 						mask = (mask + poly2mask(tmpRois{iBlood}(:,1), tmpRois{iBlood}(:,2), sz(1), sz(2))) > 0;
% 
% 					otherwise
% 						% Arbitrary number of ROIs -- add them all up
% 						cMask = cellfun(@(x) poly2mask(x(:,1), x(:,2), sz(1), sz(2)), tmpRois, 'UniformOutput', false);
% 						mask = (sum(cat(3, cMask{:}),3)) > 0;
% 				end
% 
% 				maskedImage(~mask) = nan;
% 				imgData{currentSlice}{currentTime} = maskedImage;
% 				setappdata(jAxis, 'nonMaskedImage', nonMaskedImage);
% 			end
% 
% 			setappdata(jAxis, 'imgData',        imgData);
% 			for iRoi = 1:numel(tmpRois)
% 				UpdateRoiData(jAxis, iRoi);
% 			end
% 			setappdata(jAxis, 'currentTime',    0); % Dummy value to force a reload of the image
% 		end
% 		% Need to refresh images
% 		SetTimeSlice(hFig, currentTime, currentSlice);
% 	end

	function ToggleGridfitRoi(hFig)
		% Use gridfit to fit the image over the ROIs

		if ~exist('gridfit.m', 'file')
			return
		end
		
		hAxes = getappdata(hFig, 'hAxes');

		for jAxis = reshape(hAxes(hAxes ~= 0), 1, [])
			currentSlice     = getappdata(jAxis, 'currentSlice');
			currentTime      = getappdata(jAxis, 'currentTime');
			imgData          = getappdata(jAxis, 'imgData');
			nonMaskedImgData = getappdata(jAxis, 'nonMaskedImgData');
			roi              = getappdata(jAxis, 'roi');
			tmpRois          = roi{currentSlice}{currentTime};

			if ~isempty(nonMaskedImgData)
				% Already showing masked image, so put the stored images back
				setappdata(jAxis, 'imgData',          nonMaskedImgData);
				setappdata(jAxis, 'nonMaskedImgData', []);
			else
				% Create mask based on ROIs
				nonMaskedImgData = imgData;
				for iTime = 1:numel(imgData{currentSlice})
					if isempty(imgData{currentSlice}{iTime}) || isempty(roi{currentSlice}{iTime})
						continue
					end

					tmpImg  = imgData{currentSlice}{iTime};
					tmpRois = roi{currentSlice}{iTime};
					sz      = size(tmpImg);

					% Arbitrary number of ROIs -- add them all up
					cMask = cellfun(@(x) poly2mask(x(:,1), x(:,2), sz(1), sz(2)), tmpRois, 'UniformOutput', false);
					mask = (sum(cat(3, cMask{:}),3)) > 0;

					% Apply gridfit
					[x, y] = meshgrid(1:size(tmpImg,2),1:size(tmpImg,1));

					tmpImg = gridfit(x(mask),y(mask),double(tmpImg(mask)), 1:size(tmpImg,2), 1:size(tmpImg,1), 'smoothness', 100);
% 					imgData{currentSlice}{iTime} = tmpImg;
					imgData{currentSlice}{iTime} = imgData{currentSlice}{iTime} ./ tmpImg;
				end

				setappdata(jAxis, 'nonMaskedImgData', nonMaskedImgData);
				setappdata(jAxis, 'imgData',          imgData);
			end

			for iRoi = 1:numel(tmpRois)
				UpdateRoiData(jAxis, iRoi);
			end
			setappdata(jAxis, 'currentTime',    0); % Dummy value to force a reload of the image
		end
		% Need to refresh images
		SetTimeSlice(hFig, currentTime, currentSlice);
	end

	%% ------------------------------------------------------------
% 	function StartRoiDraw(hFig, evt)
% 		% If the click occurs on an existing ROI object, abort and call that object's callback
% 		hObj  = get(hFig, 'CurrentObject');
% 		hFcn  = get(hObj, 'ButtonDownFcn');
% 		if ~isempty(hFcn)
% 			MouseDownSelect(hFig, evt);
% 			hFcn(hObj, evt);
% 			return
% 		end
% 
% 		set(hFig, 'WindowButtonDownFcn',   @InProgressRoiClick);
% 		set(hFig, 'WindowButtonMotionFcn', @InProgressRoiMove);
% 
% 		hAxis        = get(hFig,        'CurrentAxes');
% 		staticRoi    = getappdata(hFig, 'staticRoi');
% 		prefs        = getappdata(hFig, 'prefs');
% 
% 		currentPoint = get(hAxis,        'CurrentPoint');
% 		roi          = getappdata(hAxis, 'roi');
% 		roiHandles   = getappdata(hAxis, 'roiHandles');
% 		currentSlice = getappdata(hAxis, 'currentSlice');
% 		currentTime  = getappdata(hAxis, 'currentTime');
% 
% 		SelectSubplot(hFig, hAxis);
% 
% 		% Next ROI is the first empty slot (i.e. if ROI 2 of 3 is deleted, the next ROI drawn is index 2)
% 		if isempty(roi{currentSlice}{currentTime})
% 			iRoi = [];
% 		else
% 			iRoi = find(cellfun(@(x) isempty(x), roi{currentSlice}{currentTime}), 1);
% 		end
% 
% 		if isempty(iRoi)
% 			iRoi = length(roi{currentSlice}{currentTime})+1;
% 		end
% 
% 		roi{currentSlice}{currentTime}{iRoi} = repmat(currentPoint(1,1:2), [2 1]);
% 
% % This is all handled by DoStaticAndTranscribeAnd3D() now...
% % % 		% For staticRoi, propagate it to all other time frames
% % % 		if staticRoi
% % % 			for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% % % 				roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% % % 			end
% % % 		end
% % % 
% % % 		% Create plot objects if necessary
% % % 		if isempty(roiHandles)
% % % 			roiHandles = struct('line', [],  'vertices', [],  'center', [],  'label', []);
% % % 		end
% % % 
% % % 		if iRoi > length(roiHandles)
% % % 			roiHandles(iRoi) = struct('line', [],  'vertices', [],  'center', [],  'label', []);
% % % 		end
% % % 
% % % 		roiColor = prefs.roiColors(mod(iRoi-1, length(prefs.roiColors))+1);
% % % 		if isempty(roiHandles(iRoi).line)
% % % 			roiHandles(iRoi).line   = plot(0, 0, ...
% % % 		                                 'Color',      roiColor, ...
% % % 		                                 'LineStyle',  '-', ...
% % % 		                                 'Marker',     'none');
% % % 		end
% % % 
% % % 		if isempty(roiHandles(iRoi).vertices)
% % % 			roiHandles(iRoi).vertices = plot(0, 0, ...
% % % 		                                   'Color',      roiColor, ...
% % % 		                                   'LineStyle',  'none', ...
% % % 		                                   'Marker',     '.', ...
% % % 		                                   'MarkerSize', 6);
% % % 		end
% % % 
% % % 		if isempty(roiHandles(iRoi).center)
% % % 			roiHandles(iRoi).center = plot(0, 0, ...
% % % 		                                 'Color',      roiColor, ...
% % % 		                                 'LineStyle',  'none', ...
% % % 		                                 'Marker',     '+', ...
% % % 		                                 'MarkerSize', 6);
% % % 		end
% % % 
% % % 		if isempty(roiHandles(iRoi).label) && ~prefs.hideRoiLabels
% % % 			roiHandles(iRoi).label  = text(0, 0, num2str(iRoi), ...
% % % 			                               'Color',      roiColor, ...
% % % 			                               'HitTest',    'off');
% % % 		end
% % % 
% % % 		% Set positions for ROI's plot objects
% % % 		set(roiHandles(iRoi).line,     'XData', [roi{currentSlice}{currentTime}{iRoi}(:,1); roi{currentSlice}{currentTime}{iRoi}(1,1)], ...
% % % 		                               'YData', [roi{currentSlice}{currentTime}{iRoi}(:,2); roi{currentSlice}{currentTime}{iRoi}(1,2)]);
% % % 
% % % 		set(roiHandles(iRoi).vertices, 'XData', roi{currentSlice}{currentTime}{iRoi}(1,1), ...
% % % 		                               'YData', roi{currentSlice}{currentTime}{iRoi}(1,2));
% % % 
% % % 		set(roiHandles(iRoi).center,   'XData', mean(roi{currentSlice}{currentTime}{iRoi}(:,1)), ...
% % % 		                               'YData', mean(roi{currentSlice}{currentTime}{iRoi}(:,2)));
% % % 
% % % 		if ~prefs.hideRoiLabels
% % % 			set(roiHandles(iRoi).label,    'Position', [min(roi{currentSlice}{currentTime}{iRoi}(:,1))-3 min(roi{currentSlice}{currentTime}{iRoi}(:,2))-3 0]);
% % % 		end
% 
% 		% Save back to figure
% 		setappdata(hAxis, 'roi',           roi);
% 		setappdata(hAxis, 'roiHandles',    roiHandles);
% 		setappdata(hAxis, 'inProgressRoi', iRoi);
% 
% 		DoStaticAndTranscribeAnd3D(hFig, iRoi)
% 
% 		% Apply to SimpleViewer3D, if applicable
% 		hSV3D = getappdata(hFig, 'hSV3D');
% 		if ~isempty(hSV3D) && ishandle(hSV3D)
% 			hAxes = getappdata(hFig, 'hAxes');
% 			if (hAxis == hAxes(1))
% 				setappdata(hSV3D, 'roi', roi);
% 				fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 			end
% 		end
% 
% 		% Mirror this ROI on other subplots
% 		if prefs.transcribeRoi
% 			hAxes        = getappdata(hFig, 'hAxes');
% 			for jAxes = hAxes(:)'
% 				if (jAxes == 0) || (jAxes == hAxis)
% 					continue
% 				end
% 				tRoi = getappdata(jAxes, 'roi');
% 				tRoi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 
% 				% For staticRoi, propagate it to all other time frames
% 				if staticRoi
% 					for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 						tRoi{currentSlice}{iTime}{iRoi} = tRoi{currentSlice}{currentTime}{iRoi};
% 					end
% 				end
% 
% 				setappdata(jAxes, 'roi', tRoi);
% 				UpdateRoiHandles(jAxes, iRoi);
% 			end
% 		end
% 
% 		% If it's a right click (or double click), then finish the ROI immediately
% 		if any(strcmp(get(hFig, 'SelectionType'), {'alt', 'open'}))
% 			roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi}(1:end-1,:);
% 
% 			% For staticRoi, propagate it to all other time frames
% 			if staticRoi
% 				for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 					roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 				end
% 			end
% 
% 			set(roiHandles(iRoi).center,   'ButtonDownFcn', @StartRoiMove);
% 			set(roiHandles(iRoi).vertices, 'ButtonDownFcn', @VertexButtonDown);
% 			set(hFig, 'WindowButtonDownFcn', @StartRoiDraw);
% 			set(hFig, 'WindowButtonMotionFcn', '');
% 			
% 			setappdata(hAxis, 'roi',           roi);
% 			setappdata(hAxis, 'roiHandles',    roiHandles);
% 			setappdata(hAxis, 'inProgressRoi', []);
% 
% 			UpdateRoiData(hFig, iRoi);
% 
% 			% Apply to SimpleViewer3D, if applicable
% 			hSV3D = getappdata(hFig, 'hSV3D');
% 			if ~isempty(hSV3D) && ishandle(hSV3D)
% 				hAxes = getappdata(hFig, 'hAxes');
% 				if (hAxis == hAxes(1))
% 					setappdata(hSV3D, 'roi', roi);
% 					fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 					fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 				end
% 			end
% 		end
% 	end
% 
% 	function InProgressRoiClick(hFig, evt)
% 	% - Add vertex or finish ROI as requested
% 		hAxis        = get(hFig,         'CurrentAxes');
% 		currentPoint = get(hAxis,        'CurrentPoint');
% 		roi          = getappdata(hAxis, 'roi');
% 		roiHandles   = getappdata(hAxis, 'roiHandles');
% 		prefs        = getappdata(hFig,  'prefs');
% 		iRoi         = getappdata(hAxis, 'inProgressRoi');
% 		currentSlice = getappdata(hAxis, 'currentSlice');
% 		currentTime  = getappdata(hAxis, 'currentTime');
% 		staticRoi    = getappdata(hFig,  'staticRoi');
% 
% 		switch(get(hFig, 'SelectionType'))
% 			case 'normal'
% 				% Left-click: Add point
% 
% 				% Don't allow duplicate vertices
% 				if sum(abs(roi{currentSlice}{currentTime}{iRoi}(end-1,:) - currentPoint(1,1:2))) == 0
% 					return
% 				end
% 
% 				roi{currentSlice}{currentTime}{iRoi} = cat(1, roi{currentSlice}{currentTime}{iRoi}(1:end-1,:), currentPoint(1,1:2), currentPoint(1,1:2));
% 
% 				% For staticRoi, propagate it to all other time frames
% 				if staticRoi
% 					for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 						roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 					end
% 				end
% 
% 				% Add plot object for vertex if necessary
% 				iVertex = size(roi{currentSlice}{currentTime}{iRoi},1) - 1;
% 				if length(roiHandles(iRoi).vertices) < iVertex
% 					roiHandles(iRoi).vertices(iVertex) = plot(0, 0, ...
% 					                                     'Color',      prefs.roiColors(mod(iRoi-1, length(prefs.roiColors))+1), ...
% 					                                     'LineStyle',  'none', ...
% 					                                     'Marker',     '.', ...
% 					                                     'MarkerSize', 6);
% 				end
% 				set(roiHandles(iRoi).vertices(iVertex), 'XData', roi{currentSlice}{currentTime}{iRoi}(end,1), ...
% 		                                            'YData', roi{currentSlice}{currentTime}{iRoi}(end,2));
% 				setappdata(hAxis, 'roi', roi)
% 				setappdata(hAxis, 'roiHandles', roiHandles)
% 
% 				% Apply to SimpleViewer3D, if applicable
% 				hSV3D = getappdata(hFig, 'hSV3D');
% 				if ~isempty(hSV3D) && ishandle(hSV3D)
% 					hAxes = getappdata(hFig, 'hAxes');
% 					if (hAxis == hAxes(1))
% 						setappdata(hSV3D, 'roi', roi);
% 						fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 						fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 					end
% 				end
% 			case {'alt', 'open'}
% 				% Right-click or double-click: Finish ROI
% 				FinishRoi(hFig, false);
% 			case 'extend'
% 				% Shift-click to finish ROI as a line segment
% 				FinishRoi(hFig, true);
% 		end
% 	end
% 
% 	function InProgressRoiDeleteVertex(hFig, evt)
% 	% - Delete vertex from in-progress ROI
% 		hAxis        = get(hFig,         'CurrentAxes');
% 		roi          = getappdata(hAxis, 'roi');
% 		roiHandles   = getappdata(hAxis, 'roiHandles');
% 		iRoi         = getappdata(hAxis, 'inProgressRoi');
% 		currentSlice = getappdata(hAxis, 'currentSlice');
% 		currentTime  = getappdata(hAxis, 'currentTime');
% 		staticRoi    = getappdata(hFig,  'staticRoi');
% 		prefs        = getappdata(hFig,  'prefs');
% 
% 		if size(roi{currentSlice}{currentTime}{iRoi},1) < 3
% 			% --- Only a single point in the ROI -- delete the whole thing -----------
% 			if staticRoi
% 				aTimes = 1:numel(roi{currentSlice});
% 			else
% 				aTimes = currentTime;
% 			end
% 
% 			for jTime = aTimes
% 				roi{    currentSlice}{jTime}{iRoi} = [];
% 				roiMask{currentSlice}{jTime}{iRoi} = [];  % TODO: is this really necessary?
% 				roiData{currentSlice}{jTime}{iRoi} = [];  % TODO: is this really necessary?
% 			end
% 
% 			set(roiHandles(iRoi).line,     'XData', [],   'YData', []);
% 			set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
% 			set(roiHandles(iRoi).center,   'XData', [],   'YData', []);
% 			if ~prefs.hideRoiLabels
% 				set(roiHandles(iRoi).label,    'Position', [0 0 0]);
% 			end
% 
% 			setappdata(hAxis, 'inProgressRoi', []);
% 			set(hFig, 'WindowButtonDownFcn', @StartRoiDraw);
% 			set(hFig, 'WindowButtonMotionFcn', '');
% 		else
% 			% --- Delete last vertex in the ROI --------------------------------------
% 			% The last vertex is actually a placeholder for the current point
% 			roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi}([1:end-2 end],:);
% 
% 			% For staticRoi, propagate it to all other time frames
% 			if staticRoi
% 				for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 					roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 				end
% 			end
% 
% 			% Update plot objects
% 			iVertex = size(roi{currentSlice}{currentTime}{iRoi},1);
% 			set(roiHandles(iRoi).vertices(iVertex), 'XData', [],   'YData', []);
% 
% 			if ~strcmp(prefs.interpAlgorithm, 'linear')
% 				tmpInterpRoi = CubicInterp(roi{currentSlice}{currentTime}{iRoi}, prefs.interpNPts, prefs.interpAlgorithm);
% 				set(roiHandles(iRoi).line, 'XData', tmpInterpRoi(:,1), 'YData', tmpInterpRoi(:,2));
% 			else
% 				set(roiHandles(iRoi).line, 'XData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],1), ...
% 				                           'YData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],2));
% 			end
% 
% 			set(roiHandles(iRoi).center, ...
% 					'XData', mean(roi{currentSlice}{currentTime}{iRoi}(:,1)), ...
% 					'YData', mean(roi{currentSlice}{currentTime}{iRoi}(:,2)));
% 		end
% 
% 		setappdata(hAxis, 'roi', roi)
% 		setappdata(hAxis, 'roiHandles', roiHandles)
% 
% 		% Apply to SimpleViewer3D, if applicable
% 		hSV3D = getappdata(hFig, 'hSV3D');
% 		if ~isempty(hSV3D) && ishandle(hSV3D)
% 			hAxes = getappdata(hFig, 'hAxes');
% 			if (hAxis == hAxes(1))
% 				setappdata(hSV3D, 'roi', roi);
% 				fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 			end
% 		end
% 	end
% 
% 	function InProgressRoiMove(hFig, evt)
% 	% - Update plot objects for current ROI
% 		hAxis = get(hFig, 'CurrentAxes');
% 
% 		currentPoint = get(hAxis, 'CurrentPoint');
% 		roi          = getappdata(hAxis, 'roi');
% 		roiHandles   = getappdata(hAxis, 'roiHandles');
% 		iRoi         = getappdata(hAxis, 'inProgressRoi');
% 		currentSlice = getappdata(hAxis, 'currentSlice');
% 		currentTime  = getappdata(hAxis, 'currentTime');
% 		prefs        = getappdata(hFig, 'prefs');
% 
% 		% This shouldn't happen...
% 		if isempty(iRoi)
% 			return
% 		end
% 		
% 		roi{currentSlice}{currentTime}{iRoi} = cat(1, roi{currentSlice}{currentTime}{iRoi}(1:end-1,:), currentPoint(1,1:2));
% 
% 		if ~strcmp(prefs.interpAlgorithm, 'linear')
% 			tmpInterpRoi = CubicInterp(roi{currentSlice}{currentTime}{iRoi}, prefs.interpNPts, prefs.interpAlgorithm);
% 			set(roiHandles(iRoi).line, 'XData', tmpInterpRoi(:,1), 'YData', tmpInterpRoi(:,2));
% 		else
% 			set(roiHandles(iRoi).line, 'XData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],1), ...
% 			                           'YData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],2));
% 		end
% 
% 		set(roiHandles(iRoi).center, ...
% 		    'XData', mean(roi{currentSlice}{currentTime}{iRoi}(:,1)), ...
% 		    'YData', mean(roi{currentSlice}{currentTime}{iRoi}(:,2)));
% 
% 		if ~prefs.hideRoiLabels
% 			set(roiHandles(iRoi).label, ...
% 			    'Position', [min(roi{currentSlice}{currentTime}{iRoi}(:,1))-3 min(roi{currentSlice}{currentTime}{iRoi}(:,2))-3 0]);
% 		end
% 
% 		% Vertex drawing is deliberately omitted to provide better visualization of the to-be-placed vertex
% 		setappdata(hAxis, 'roi', roi);
% 
% 		% Mirror this ROI on other subplots
% 		if prefs.transcribeRoi
% 			hAxes = getappdata(hFig, 'hAxes');
% 			for jAxes = hAxes(:)'
% 				if (jAxes == 0) || (jAxes == hAxis)
% 					continue
% 				end
% 				tRoi = getappdata(jAxes, 'roi');
% 				tRoi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 				setappdata(jAxes, 'roi', tRoi);
% 				UpdateRoiHandles(jAxes, iRoi);
% 			end
% 		end
% 
% 		% Apply to SimpleViewer3D, if applicable
% 		hSV3D = getappdata(hFig, 'hSV3D');
% 		if ~isempty(hSV3D) && ishandle(hSV3D)
% 			hAxes = getappdata(hFig, 'hAxes');
% 			if (hAxis == hAxes(1))
% 				setappdata(hSV3D, 'roi', roi);
% 				fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 			end
% 		end
% 	end
% 
% 	function FinishRoi(hFig, bIsLine)
% 	% - Draw last vertex, update ROI data and 'inProgressRoi'
% 	% - Set the mouse callbacks appropriately
% 		hAxis        = get(hFig,         'CurrentAxes');
% 		currentPoint = get(hAxis,        'CurrentPoint');
% 		roi          = getappdata(hAxis, 'roi');
% 		roiHandles   = getappdata(hAxis, 'roiHandles');
% 		iRoi         = getappdata(hAxis, 'inProgressRoi');
% 		prefs        = getappdata(hFig,  'prefs');
% 		currentSlice = getappdata(hAxis, 'currentSlice');
% 		currentTime  = getappdata(hAxis, 'currentTime');
% 		staticRoi    = getappdata(hFig,  'staticRoi');
% 
% 		if bIsLine
% 			roi{currentSlice}{currentTime}{iRoi} = [roi{currentSlice}{currentTime}{iRoi}(1:end-1,:); currentPoint(1,1:2); nan(1,2)];
% 		else
% 			roi{currentSlice}{currentTime}{iRoi} = [roi{currentSlice}{currentTime}{iRoi}(1:end-1,:); currentPoint(1,1:2)];
% 		end
% 
% 		% For staticRoi, propagate it to all other time frames
% 		if staticRoi
% 			for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 				roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 			end
% 		end
% 
% 		% Add plot object for vertex if necessary
% % 		iVertex = size(roi{currentSlice}{currentTime}{iRoi},1);
% % 		if length(roiHandles(iRoi).vertices) < iVertex
% % 			roiHandles(iRoi).vertices(iVertex) = plot(0, 0, ...
% % 			                                     'Color',      prefs.roiColors(mod(iRoi-1, length(prefs.roiColors))+1), ...
% % 			                                     'LineStyle',  'none', ...
% % 			                                     'Marker',     '.', ...
% % 			                                     'MarkerSize', 6);
% % 		end
% % 		set(roiHandles(iRoi).vertices(iVertex), 'XData', roi{currentSlice}{currentTime}{iRoi}(end,1), ...
% % 		                                        'YData', roi{currentSlice}{currentTime}{iRoi}(end,2));
% 
% 		% New code that allows more than 1 point to be added here, as is necessary
% 		% for lines
% 		for iVertex = 1:size(roi{currentSlice}{currentTime}{iRoi},1)
% 			if iVertex > length(roiHandles(iRoi).vertices)
% 				roiHandles(iRoi).vertices(iVertex) = plot(0, 0, ...
% 				                                     'Color',      prefs.roiColors(mod(iRoi-1, length(prefs.roiColors))+1), ...
% 				                                     'LineStyle',  'none', ...
% 				                                     'Marker',     '.', ...
% 				                                     'MarkerSize', 6);
% 			end
% 		end
% 
% 		set(roiHandles(iRoi).vertices(iVertex), 'XData', roi{currentSlice}{currentTime}{iRoi}(end,1), ...
% 		                                        'YData', roi{currentSlice}{currentTime}{iRoi}(end,2));
% 
% 		% Need to do end-1 for lines, as end is simply a placeholder nan
% 		if bIsLine
% 			set(roiHandles(iRoi).vertices(iVertex), 'XData', roi{currentSlice}{currentTime}{iRoi}(end-1,1), ...
% 			                                        'YData', roi{currentSlice}{currentTime}{iRoi}(end-1,2));
% 		end
% 
% 
% 
% 		setappdata(hAxis, 'roi',           roi);
% 		setappdata(hAxis, 'roiHandles',    roiHandles);
% 		setappdata(hAxis, 'inProgressRoi', []);
% 
% 		set(roiHandles(iRoi).center,   'ButtonDownFcn', @StartRoiMove);
% 		set(roiHandles(iRoi).vertices, 'ButtonDownFcn', @VertexButtonDown);
% 		set(roiHandles(iRoi).line,     'ButtonDownFcn', @AddVertexFromEdge);
% 
% 		set(hFig, 'WindowButtonDownFcn', @StartRoiDraw);
% 		set(hFig, 'WindowButtonMotionFcn', '');
% 
% 		UpdateRoiData(hFig, iRoi);
% 
% 		% Automatically interpolate vertices
% 		if prefs.autoInterpolateRois
% 			InterpolateRoi(hAxis, iRoi);
% 		end
% 
% 		% Apply to SimpleViewer3D, if applicable
% 		hSV3D = getappdata(hFig, 'hSV3D');
% 		if ~isempty(hSV3D) && ishandle(hSV3D)
% 			hAxes = getappdata(hFig, 'hAxes');
% 			if (hAxis == hAxes(1))
% 				setappdata(hSV3D, 'roi', roi);
% 				fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 			end
% 		end
% 	end

	function StartRoiDraw(hFig, evt)
		% If the click occurs on an existing ROI object, abort and call that object's callback
		hObj  = get(hFig, 'CurrentObject');
		hFcn  = get(hObj, 'ButtonDownFcn');
		if ~isempty(hFcn)
			MouseDownSelect(hFig, evt);
			hFcn(hObj, evt);
			return
		end

		CurrentAxes  = get(hFig,               'CurrentAxes');		
		% Make sure the CurrentAxes corresponds to a subplot
		hAxes = getappdata(hFig, 'hAxes');
		if ~any(hAxes == CurrentAxes)
			return
		end

		set(hFig, 'WindowButtonDownFcn',   @InProgressRoiClick);
		set(hFig, 'WindowButtonMotionFcn', @InProgressRoiMove);
		
		currentPoint = get(CurrentAxes,        'CurrentPoint');
		roi          = getappdata(CurrentAxes, 'roi');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');

		SelectSubplot(hFig, CurrentAxes);

		% Next ROI is the first empty slot (i.e. if ROI 2 of 3 is deleted, the next ROI drawn is index 2)
		if isempty(roi{currentSlice}{currentTime})
			iRoi = [];
		else
			iRoi = find(cellfun(@(x) isempty(x), roi{currentSlice}{currentTime}), 1);
		end

		if isempty(iRoi)
			iRoi = length(roi{currentSlice}{currentTime})+1;
		end
		setappdata(CurrentAxes, 'inProgressRoi',     iRoi);
		setappdata(hFig,        'inProgressRoiAxis', CurrentAxes);

		% Add the data point, plus a 'repmat'ed point to follow the cursor for preview
		roi{currentSlice}{currentTime}{iRoi} = repmat(currentPoint(1,1:2), [2 1]);
		bDontUpdateData = true;

		% If it's a right click (or double click), then finish the ROI immediately
		if any(strcmp(get(hFig, 'SelectionType'), {'alt', 'open'}))
			roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi}(1:end-1,:);
			set(hFig, 'WindowButtonDownFcn',   @StartRoiDraw);
			set(hFig, 'WindowButtonMotionFcn', '');
			setappdata(CurrentAxes, 'inProgressRoi',     []);
			setappdata(hFig,        'inProgressRoiAxis', []);
			bDontUpdateData = false;
		end
		setappdata(CurrentAxes, 'roi', roi);

		DoStaticAndTranscribeAnd3D(hFig, iRoi, bDontUpdateData);

% 		% Hide cursor while drawing to make things easier to see
% 		set(hFig, 'Pointer', 'custom', 'PointerShapeCData', nan(16))
	end

	%% ------------------------------------------------------------
	function InProgressRoiClick(hFig, evt)
	% - Add vertex or finish ROI as requested
		CurrentAxes  = get(hFig,               'CurrentAxes');
		currentPoint = get(CurrentAxes,        'CurrentPoint');
		roi          = getappdata(CurrentAxes, 'roi');
		iRoi         = getappdata(CurrentAxes, 'inProgressRoi');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');

		% Abort if the subplot changed
		inProgressRoiAxis = getappdata(hFig, 'inProgressRoiAxis');
		if (CurrentAxes ~= inProgressRoiAxis)
			% Actually, change it back
			set(hFig, 'CurrentAxes', inProgressRoiAxis);
			return
		end

		switch(get(hFig, 'SelectionType'))
			case 'normal'
				% Left-click: Add point

				% Don't allow duplicate vertices
				if sum(abs(roi{currentSlice}{currentTime}{iRoi}(end-1,:) - currentPoint(1,1:2))) == 0
					return
				end

				roi{currentSlice}{currentTime}{iRoi} = cat(1, roi{currentSlice}{currentTime}{iRoi}(1:end-1,:), currentPoint(1,1:2), currentPoint(1,1:2));
				setappdata(CurrentAxes, 'roi', roi);
				DoStaticAndTranscribeAnd3D(hFig, iRoi, true);
			case {'alt'}%, 'open'}
				% Right-click or double-click: Finish ROI
				FinishRoi(hFig, false);
			case 'extend'
				% Shift-click to finish ROI as a line segment
				FinishRoi(hFig, true);
			case 'open'
				% Double-click to finish as an ROI without interpolating
				FinishRoi(hFig, false, true);
		end
	end

	%% ------------------------------------------------------------
	function InProgressRoiDeleteVertex(hFig, evt)
	% - Delete vertex from in-progress ROI
		CurrentAxes  = get(hFig,               'CurrentAxes');
		roi          = getappdata(CurrentAxes, 'roi');
		roiHandles   = getappdata(CurrentAxes, 'roiHandles');
		iRoi         = getappdata(CurrentAxes, 'inProgressRoi');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');

		if size(roi{currentSlice}{currentTime}{iRoi},1) < 3
			% Only a single point in the ROI -- delete the whole thing
			roi{currentSlice}{currentTime}{iRoi} = [];

			setappdata(CurrentAxes, 'inProgressRoi',     []);
			setappdata(hFig,        'inProgressRoiAxis', []);
			set(hFig, 'WindowButtonDownFcn',   @StartRoiDraw);
			set(hFig, 'WindowButtonMotionFcn', '');
		else
			% Delete last vertex in the ROI
			% The last vertex is actually a placeholder for the current point
			roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi}([1:end-2 end],:);
		end
		setappdata(CurrentAxes, 'roi', roi)
		DoStaticAndTranscribeAnd3D(hFig, iRoi, true);
		
		% If we're still drawing, hide the vertex of the preview ROI on the last point
		% Consider doing this for other subplots too
		if size(roi{currentSlice}{currentTime}{iRoi},1) >= 3
			iVertex = size(roi{currentSlice}{currentTime}{iRoi},1);
			set(roiHandles(iRoi).vertices(iVertex), 'XData', [])
			set(roiHandles(iRoi).vertices(iVertex), 'YData', [])
		end
	end

	%% ------------------------------------------------------------
	function InProgressRoiMove(hFig, evt)
	% - Update plot objects for current ROI
		CurrentAxes  = getappdata(hFig,        'inProgressRoiAxis');
% 		CurrentAxes  = get(hFig,               'CurrentAxes');
		currentPoint = get(CurrentAxes,        'CurrentPoint');
		roi          = getappdata(CurrentAxes, 'roi');
		roiHandles   = getappdata(CurrentAxes, 'roiHandles');
		currRoi      = getappdata(CurrentAxes, 'inProgressRoi');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');

		% This shouldn't happen...
		if isempty(currRoi)
			return
		end

		roi{currentSlice}{currentTime}{currRoi} = cat(1, roi{currentSlice}{currentTime}{currRoi}(1:end-1,:), currentPoint(1,1:2));
		setappdata(CurrentAxes, 'roi', roi)
		
		% This may be sub-optimal, as it re-draws everything.  Consider an optimized version
% 		DoStaticAndTranscribeAnd3D(hFig, iRoi, true);

		% --------------------------------------------------------------------------
		% This is an version of DoStaticAndTranscribeAnd3D, optimized for speed with
		% the knowledge that only the last point is moving		
		staticRoi       = getappdata(hFig, 'staticRoi');
		prefs           = getappdata(hFig, 'prefs');
		hAxes           = getappdata(hFig, 'hAxes');

		% Figure out which axes are getting updated
		if iscell(prefs.transcribeRoi)
			currentSubplot = find(hAxes == CurrentAxes);
			goodInds = cellfun(@(x) any(x == currentSubplot), prefs.transcribeRoi);
			if any(goodInds)
				aSubplots = prefs.transcribeRoi{goodInds};
				aAxes     = hAxes(aSubplots);
			else
				aAxes = CurrentAxes;
			end
		elseif prefs.transcribeRoi
			aAxes = hAxes(hAxes ~= 0);
		elseif ~prefs.transcribeRoi
			aAxes = CurrentAxes;
		end
		aAxes = reshape(aAxes, 1, []);

		% Consider not doing this if this takes a while.  The actual ROI doesn't
		% necessary need to be updated until the ROI is finished
		for jAxis = aAxes
			tRoi = getappdata(jAxis, 'roi');

			% Already did data for current one
			if (jAxis ~= CurrentAxes)
				% Range check (these are clipped, not cyclic)
				tSlice = min([currentSlice length(tRoi        )]);
				tTime  = min([currentTime  length(tRoi{tSlice})]);

				tRoi{tSlice}{tTime}(currRoi) = roi{currentSlice}{currentTime}(currRoi);

				% For staticRoi, propagate it to all other time frames
				if staticRoi
					for iTime = [1:tTime-1 tTime+1:numel(tRoi{tSlice})]
						tRoi{tSlice}{iTime}(currRoi) = tRoi{tSlice}{tTime}(currRoi);
					end
				end
				setappdata(jAxis, 'roi', tRoi);
			else
				tSlice = currentSlice;
				tTime  = currentTime;
			end
			% --------------------------------------------------------------------------
			% This is an version of UpdateRoiHandles, optimized for speed with
			% the knowledge that only the last point is moving		
			roiHandles   = getappdata(jAxis, 'roiHandles');

			% Update plot object positions
			if ~strcmp(prefs.interpAlgorithm, 'linear') && any(strcmp(prefs.interpRois, {'display', 'finish', 'always'}))
				tmpInterpRoi = CubicInterp(tRoi{tSlice}{tTime}{currRoi}, [], prefs.interpAlgorithm);
				set(roiHandles(currRoi).line, 'XData', tmpInterpRoi(:,1), 'YData', tmpInterpRoi(:,2));
			else
				% Line segments must be plotted differently as of R2014b
				if verLessThan('matlab','8.4.0')
					set(roiHandles(currRoi).line, 'XData', tRoi{tSlice}{tTime}{currRoi}([1:end 1],1), ...
					                              'YData', tRoi{tSlice}{tTime}{currRoi}([1:end 1],2));
				else
					roiHandles(currRoi).line.XData = tRoi{tSlice}{tTime}{currRoi}([1:end 1],1);
					roiHandles(currRoi).line.YData = tRoi{tSlice}{tTime}{currRoi}([1:end 1],2);
				end
			end

			set(roiHandles(currRoi).center, 'XData',    mean(tRoi{tSlice}{tTime}{currRoi}(:,1)), ...
			                                'YData',    mean(tRoi{tSlice}{tTime}{currRoi}(:,2)));

			% Don't bother updating the last vertex for the current axis, as we were going to hide it anyways
			if (jAxis ~= CurrentAxes)

				for iVertex = size(tRoi{tSlice}{tTime}{currRoi},1) % Just do the last one, not all. i.e. 1:size(tRoi{currentSlice}{currentTime}{currRoi},1)
	% % 				if verLessThan('matlab','8.4.0')
	% % 					roiHandles(iRoi).vertices(iVertex).XData = tRoi{currentSlice}{currentTime}{iRoi}(iVertex,1);
	% % 					roiHandles(iRoi).vertices(iVertex).YData = tRoi{currentSlice}{currentTime}{iRoi}(iVertex,2);
	% % 				else
						set(roiHandles(currRoi).vertices(iVertex), 'XData', tRoi{tSlice}{tTime}{currRoi}(iVertex,1), 'YData', tRoi{tSlice}{tTime}{currRoi}(iVertex,2));
	% % 				end
				end
			end
		end

		% Apply to SimpleViewer3D, if applicable
		hSV3D = getappdata(hFig, 'hSV3D');
		if ~isempty(hSV3D) && ishandle(hSV3D)
			if (CurrentAxes == hAxes(1))
				setappdata(hSV3D, 'roi', getappdata(CurrentAxes, 'roi'));
				fcnHandles = getappdata(hSV3D, 'fcnHandles');
				% FIXME: These functions should really support multiple selectedRois...
				for iRoi = selectedRois
					fcnHandles.UpdateRoiHandles(hSV3D, currRoi);
				end
			end
		end

% 		% Hide the last vertex to provide better visualization of the to-be-placed vertex
% 		% Consider doing this for other subplots too
% 		iVertex = size(roi{currentSlice}{currentTime}{currRoi},1);
% 		set(roiHandles(currRoi).vertices(iVertex), 'XData', [])
% 		set(roiHandles(currRoi).vertices(iVertex), 'YData', [])
	end

	%% ------------------------------------------------------------
	function FinishRoi(hFig, bIsLine, bSkipInterp)
	% - Draw last vertex, update ROI data and 'inProgressRoi'
	% - Set the mouse callbacks appropriately
		CurrentAxes  = get(hFig,               'CurrentAxes');
		currentPoint = get(CurrentAxes,        'CurrentPoint');
		roi          = getappdata(CurrentAxes, 'roi');
		roiHandles   = getappdata(CurrentAxes, 'roiHandles');
		iRoi         = getappdata(CurrentAxes, 'inProgressRoi');
		prefs        = getappdata(hFig,        'prefs');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');
		currentTime  = getappdata(CurrentAxes, 'currentTime');
		staticRoi    = getappdata(hFig,        'staticRoi');

		tol = 1e-6; % distances smaller than this are considered negligble
		
		if ~exist('bSkipInterp', 'var') || isempty(bSkipInterp)
			bSkipInterp = false;
		end

		if bIsLine
			roi{currentSlice}{currentTime}{iRoi} = [roi{currentSlice}{currentTime}{iRoi}(1:end-1,:); currentPoint(1,1:2); nan(1,2)];
		else
			roi{currentSlice}{currentTime}{iRoi} = [roi{currentSlice}{currentTime}{iRoi}(1:end-1,:); currentPoint(1,1:2)];
		end
		
		% Remove duplicate end point that arises from a double-click finish
		if (sum(abs(diff(roi{currentSlice}{currentTime}{iRoi}(end-1:end,:)))) < tol)
			roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi}([1:end-2 end],:);
		end

		setappdata(CurrentAxes, 'roi',           roi);

% 		% Need to do end-1 for lines, as end is simply a placeholder nan
% 		if bIsLine
% 			set(roiHandles(iRoi).vertices(iVertex), 'XData', roi{currentSlice}{currentTime}{iRoi}(end-1,1), ...
% 			                                        'YData', roi{currentSlice}{currentTime}{iRoi}(end-1,2));
% 		end

		setappdata(CurrentAxes, 'roiHandles',        roiHandles);
		setappdata(CurrentAxes, 'inProgressRoi',     []);
		setappdata(hFig,        'inProgressRoiAxis', []);
		set(hFig, 'WindowButtonDownFcn', @StartRoiDraw);
		set(hFig, 'WindowButtonMotionFcn', '');

		if ~bSkipInterp && (prefs.autoInterpolateRois || any(strcmp(prefs.interpRois, {'finish', 'always'})))
			% Automatically interpolate vertices (this will call DoStaticAndTranscribeAnd3D)
			InterpolateRoi(CurrentAxes, iRoi);
		else
			DoStaticAndTranscribeAnd3D(hFig, iRoi);
		end
		
		% Update ROI data
		UpdateRoiData(hFig, iRoi, 1:numel(roi{currentSlice}))
		
% 		% Un-hide cursor
% 		setptr(hFig, 'datacursor');
	end

  % ------------------------------------------------------------
	function AddVertexFromEdge(src, evt)
		% Don't process this if it's passed on from another subfunction (otherwise
		% it'll be called twice)
		if numel(dbstack) > 1
			return
		end
		
		% TODO: this suggests that we could replace individual vertices with just
		% the line and using IntersectionPoint to determine the closest vertex
% 		evt.IntersectionPoint
		
		hFig  = get(get(src, 'Parent'), 'Parent');

		% Shift-click to add vertices
		if ~strcmp(get(hFig, 'SelectionType'), 'extend')
			return
		end

		% Don't allow this when nudging, as it gets in the way and we already have a way of adjusting the ROI
		if any(strcmp(getappdata(hFig, 'strMode'), {'Nudge', 'Circle'}))
			return
		end

		% Also not allowed while drawing
		if ~isempty(getappdata(get(hFig, 'CurrentAxes'), 'inProgressRoi'))
			return
		end

		hAxis        = get(hFig,         'CurrentAxes');
		currentPoint = get(hAxis,        'CurrentPoint');
		roi          = getappdata(hAxis, 'roi');
		roiHandles   = getappdata(hAxis, 'roiHandles');
		selectedRoi  = getappdata(hFig,  'selectedRoi');
		prefs        = getappdata(hFig,  'prefs');
		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');
		staticRoi    = getappdata(hFig,  'staticRoi');

		if isempty(selectedRoi)
			return
		end

		% Make sure the line belongs to the selectedRoi.  (Odd case can happen when
		% finishing another ROI with a shift-click)
		if (src ~= roiHandles(selectedRoi).line)
			return
		end
		
		% Determine after which vertex the new one goes:
		% Find the area of the new vertex with each pair of adjacent vertices
		% The minimal area set is the most co-linear, and thus the vertex goes here
		% Determinant is proportional to area (otherwise use polyarea, which is slower)
		% Note: This can fail when the new vertex should lie between A and B, but
		% the distance between A and B is much greater than the distance between B
		% and C.  Therefore, normalize to the distance of the two test vertices.
		currentRoiVertices = [roi{currentSlice}{currentTime}{selectedRoi}; roi{currentSlice}{currentTime}{selectedRoi}(1,:)];

		aNormDet = realmax*ones(1,size(currentRoiVertices,1)-1);
		for iInd = 1:size(currentRoiVertices,1)-1
			XData = [currentRoiVertices(iInd:iInd+1,1); currentPoint(1,1)];
			YData = [currentRoiVertices(iInd:iInd+1,2); currentPoint(1,2)];

			dist = sqrt(diff(XData(1:2))^2 + diff(YData(1:2))^2);
			tDet = abs(det([XData(1)-XData(3)  YData(1)-YData(3); ...
			                XData(2)-XData(3)  YData(2)-YData(3)]));
			
			aNormDet(iInd) = tDet/dist;
			
		end
		newInd = find(aNormDet == min(aNormDet));

		% Add vertex to ROI
		roi{currentSlice}{currentTime}{selectedRoi} = cat(1,roi{currentSlice}{currentTime}{selectedRoi}(1:newInd,:), currentPoint(1,1:2), roi{currentSlice}{currentTime}{selectedRoi}(newInd+1:end,:));

		% For staticRoi, propagate it to all other time frames
		if staticRoi
			for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
				roi{currentSlice}{iTime}{selectedRoi} = roi{currentSlice}{currentTime}{selectedRoi};
			end
		end

		% Update the line
		if ~strcmp(prefs.interpAlgorithm, 'linear') && any(strcmp(prefs.interpRois, {'display', 'always'}))
			tmpInterpRoi = CubicInterp(roi{currentSlice}{currentTime}{selectedRoi}, [], prefs.interpAlgorithm);
			set(roiHandles(selectedRoi).line, 'XData', tmpInterpRoi(:,1), 'YData', tmpInterpRoi(:,2));
		else
			set(roiHandles(selectedRoi).line, 'XData', roi{currentSlice}{currentTime}{selectedRoi}([1:end 1],1), ...
			                                  'YData', roi{currentSlice}{currentTime}{selectedRoi}([1:end 1],2));
		end

		% Update the ceter
		set(roiHandles(selectedRoi).center, ...
		    'XData', mean(roi{currentSlice}{currentTime}{selectedRoi}(:,1)), ...
		    'YData', mean(roi{currentSlice}{currentTime}{selectedRoi}(:,2)));

		% Draw the vertex
		hNewVertex = plot(currentPoint(1,1), currentPoint(1,2), ...
					            'Color',      prefs.roiColors(mod(selectedRoi-1, size(prefs.roiColors,1))+1,:), ...
					            'LineStyle',  'none', ...
					            'Marker',     '+');%, ...
% 					            'MarkerSize', 6);
		set(hNewVertex, 'ButtonDownFcn', @VertexButtonDown);
		roiHandles(selectedRoi).vertices = cat(2, roiHandles(selectedRoi).vertices(1:newInd), hNewVertex, roiHandles(selectedRoi).vertices(newInd+1:end));


		setappdata(hAxis, 'roi', roi)
		setappdata(hAxis, 'roiHandles', roiHandles)

		% Apply to SimpleViewer3D, if applicable
		hSV3D = getappdata(hFig, 'hSV3D');
		if ~isempty(hSV3D) && ishandle(hSV3D)
			hAxes = getappdata(hFig, 'hAxes');
			if (hAxis == hAxes(1))
				setappdata(hSV3D, 'roi', roi);
				fcnHandles = getappdata(hSV3D, 'fcnHandles');
				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
			end
		end
	end

	%% ------------------------------------------------------------
	function DeleteVertex(hFig, hObj)
	% - Delete vertex from a finished ROI
		hAxis        = get(hFig,         'CurrentAxes');
		roi          = getappdata(hAxis, 'roi');
		roiHandles   = getappdata(hAxis, 'roiHandles');
		selectedRoi  = getappdata(hFig,  'selectedRoi');
		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');
		staticRoi    = getappdata(hFig,  'staticRoi');
		prefs        = getappdata(hFig,  'prefs');

		% Find the clicked vertex
		iVertex = find(roiHandles(selectedRoi).vertices == hObj);

		% If it's the only vertex, delete the whole thing
		% Note, however, that this never happens because the center + is on top and so the vertex is unclickable
		if size(roi{currentSlice}{currentTime}{selectedRoi},1) == 1
			DeleteRois(hFig, getappdata(hFig, 'selectedRoi'))
			return
		end

		% Remove the vertex
		roi{currentSlice}{currentTime}{selectedRoi} = roi{currentSlice}{currentTime}{selectedRoi}([1:iVertex-1 iVertex+1:end],:);

		% For staticRoi, propagate it to all other time frames
		if staticRoi
			for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
				roi{currentSlice}{iTime}{selectedRoi} = roi{currentSlice}{currentTime}{selectedRoi};
			end
		end

		% Update the line
		if ~strcmp(prefs.interpAlgorithm, 'linear') && any(strcmp(prefs.interpRois, {'display', 'always'}))
			tmpInterpRoi = CubicInterp(roi{currentSlice}{currentTime}{selectedRoi}, [], prefs.interpAlgorithm);
			set(roiHandles(selectedRoi).line, 'XData', tmpInterpRoi(:,1), 'YData', tmpInterpRoi(:,2));
		else
			set(roiHandles(selectedRoi).line, 'XData', roi{currentSlice}{currentTime}{selectedRoi}([1:end 1],1), ...
			                                  'YData', roi{currentSlice}{currentTime}{selectedRoi}([1:end 1],2));
		end

		% Update the center
		set(roiHandles(selectedRoi).center, ...
		    'XData', mean(roi{currentSlice}{currentTime}{selectedRoi}(:,1)), ...
		    'YData', mean(roi{currentSlice}{currentTime}{selectedRoi}(:,2)));

		% Hide the plot object
		set(roiHandles(selectedRoi).vertices(iVertex), 'XData', [],   'YData', []);
		roiHandles(selectedRoi).vertices = roiHandles(selectedRoi).vertices([1:iVertex-1 iVertex+1:end iVertex]);
		
		setappdata(hAxis, 'roi', roi)
		setappdata(hAxis, 'roiHandles', roiHandles)

		% Apply to SimpleViewer3D, if applicable
		hSV3D = getappdata(hFig, 'hSV3D');
		if ~isempty(hSV3D) && ishandle(hSV3D)
			hAxes = getappdata(hFig, 'hAxes');
			if (hAxis == hAxes(1))
				setappdata(hSV3D, 'roi', roi);
				fcnHandles = getappdata(hSV3D, 'fcnHandles');
				fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
			end
		end
	end

  % ------------------------------------------------------------
	function VertexButtonDown(src, evt)
		if strcmp(get(src, 'Type'), 'figure')
			hFig = src;
		else
			hFig = get(get(src, 'Parent'), 'Parent');
		end

		% Don't allow this when nudging, as it gets in the way and we already have a way of adjusting the ROI
		if any(strcmp(getappdata(hFig, 'strMode'), {'Nudge', 'Circle'}))
			return
		end

		% Also not allowed while drawing
		if ~isempty(getappdata(get(hFig, 'CurrentAxes'), 'inProgressRoi'))
			return
		end

		if strcmp(get(hFig, 'SelectionType'), 'extend')
			% Shift-click to delete vertex
			DeleteVertex(hFig, src);
		elseif strcmp(get(hFig, 'SelectionType'), 'open')
			% Ignore double clicks, as they won't have a corresponding 'ButtonUp' to
			% know when to stop moving the vertex
		elseif strcmp(get(hFig, 'SelectionType'), 'alt')
			% Ingore right clicks too, to avoid calls when finishing an ROI
		else
			% Otherwise just move the ROI
			set(hFig, 'WindowButtonUpFcn',     @EndRoiMove);
			set(hFig, 'WindowButtonMotionFcn', @RoiMove);
		end
	end

  % ------------------------------------------------------------
	function StartRoiMove(src, evt)
		if strcmp(get(src, 'Type'), 'figure')
			hFig = src;
		else
			hFig = get(get(src, 'Parent'), 'Parent');
		end

		% Don't allow this when nudging, as it gets in the way and we already have a way of adjusting the ROI
		if any(strcmp(getappdata(hFig, 'strMode'), {'Nudge', 'Circle'}))
			return
		end

		% Also not allowed while drawing
		if ~isempty(getappdata(get(hFig, 'CurrentAxes'), 'inProgressRoi'))
			return
		end

		% Disable usage for right clicks, such as when finishing an ROI
		if strcmp(get(hFig, 'SelectionType'), 'alt')
			return
		end
		
		set(hFig, 'WindowButtonMotionFcn', @RoiMove);
		set(hFig, 'WindowButtonUpFcn',     @EndRoiMove);
	end

  % ------------------------------------------------------------
	function EndRoiMove(src, evt)
% 		hFig = get(get(src, 'Parent'), 'Parent');
		hFig = src;

		set(hFig, 'WindowButtonMotionFcn', '');
		set(hFig, 'WindowButtonUpFcn',     '');
	end

  % ------------------------------------------------------------
% 	function RoiMove(hFig, evt)
% 	% - Find the selected ROI, update the ROI data, then call UpdateRoiData on it
% 		hAxis = get(hFig, 'CurrentAxes');
% 		hObj  = get(hFig, 'CurrentObject');
% 
% 		currentPoint = get(hAxis,        'CurrentPoint');
% 		roi          = getappdata(hAxis, 'roi');
% 		roiHandles   = getappdata(hAxis, 'roiHandles');
% 		currentSlice = getappdata(hAxis, 'currentSlice');
% 		currentTime  = getappdata(hAxis, 'currentTime');
% 		staticRoi    = getappdata(hFig,   'staticRoi');
% 
% 		for iRoi = 1:length(roiHandles)
% 			% Selected object is an ROI center
% 			if (roiHandles(iRoi).center == hObj)
% 				centerDiff = currentPoint(1,1:2) - mean(roi{currentSlice}{currentTime}{iRoi},1);
% 				roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi} + repmat(centerDiff, [size(roi{currentSlice}{currentTime}{iRoi},1) 1]);
% 
% 				% For staticRoi, propagate it to all other time frames
% 				if staticRoi
% 					for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 						roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 					end
% 				end
% 
% 				setappdata(hAxis, 'roi', roi);
% 				UpdateRoiData(hFig, iRoi);
% 
% 				% Apply to SimpleViewer3D, if applicable
% 				hSV3D = getappdata(hFig, 'hSV3D');
% 				if ~isempty(hSV3D) && ishandle(hSV3D)
% 					hAxes = getappdata(hFig, 'hAxes');
% 					if (hAxis == hAxes(1))
% 						setappdata(hSV3D, 'roi', roi);
% 						fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 						fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 					end
% 				end
% 
% 				return
% 			end
% 			
% 			% Selected object is an ROI vertex
% 			iVertex = find(roiHandles(iRoi).vertices == hObj);
% 			if ~isempty(iVertex)
% 				roi{currentSlice}{currentTime}{iRoi}(iVertex,:) = currentPoint(1,1:2);
% 
% 				% For staticRoi, propagate it to all other time frames
% 				if staticRoi
% 					for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 						roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 					end
% 				end
% 
% 				setappdata(hAxis, 'roi', roi);
% 				UpdateRoiData(hFig, iRoi);
% 
% 				% Apply to SimpleViewer3D, if applicable
% 				hSV3D = getappdata(hFig, 'hSV3D');
% 				if ~isempty(hSV3D) && ishandle(hSV3D)
% 					hAxes = getappdata(hFig, 'hAxes');
% 					if (hAxis == hAxes(1))
% 						setappdata(hSV3D, 'roi', roi);
% 						fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 						fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 					end
% 				end
% 
% 				return
% 			end
% 		end
% 	end
  % ------------------------------------------------------------
	function RoiMove(hFig, evt)
	% - Find the selected ROI, update the ROI data, then call UpdateRoiData on it
		hAxis = get(hFig, 'CurrentAxes');
		hObj  = get(hFig, 'CurrentObject');

		if isempty(hObj)
			return
		end

		currentPoint = get(hAxis,        'CurrentPoint');
		roi          = getappdata(hAxis, 'roi');
		roiHandles   = getappdata(hAxis, 'roiHandles');
		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');

		for iRoi = 1:length(roiHandles)
			% Selected object is an ROI center
			if (roiHandles(iRoi).center == hObj)
				centerDiff = currentPoint(1,1:2) - mean(roi{currentSlice}{currentTime}{iRoi},1);
				roi{currentSlice}{currentTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi} + repmat(centerDiff, [size(roi{currentSlice}{currentTime}{iRoi},1) 1]);
				setappdata(hAxis, 'roi', roi);
				DoStaticAndTranscribeAnd3D(hFig, iRoi);
			end
			
			% Selected object is an ROI vertex
			iVertex = find(roiHandles(iRoi).vertices == hObj);
			if ~isempty(iVertex)
				roi{currentSlice}{currentTime}{iRoi}(iVertex,:) = currentPoint(1,1:2);
				setappdata(hAxis, 'roi', roi);
				DoStaticAndTranscribeAnd3D(hFig, iRoi);
			end
		end
	end

  % ------------------------------------------------------------
	function UpdateAllRoiHandles(hFig)
	% - Update all visible ROIs on all subplots
		hAxes = getappdata(hFig, 'hAxes');
		for hAxis = reshape(hAxes(hAxes ~= 0), 1, [])
			roi          = getappdata(hAxis, 'roi');
			currentSlice = getappdata(hAxis, 'currentSlice');
			currentTime  = getappdata(hAxis, 'currentTime');
			
			for iRoi = 1:numel(roi{currentSlice}{currentTime})
				UpdateRoiHandles(hAxis, iRoi)
			end
		end
	end

  % ------------------------------------------------------------
	function UpdateRoiHandles(hAxis, iRoi)
	% - Update mask, data, and plots for a given ROI
		roi          = getappdata(hAxis, 'roi');
		roiHandles   = getappdata(hAxis, 'roiHandles');
		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');
		prefs        = getappdata(get(hAxis, 'Parent'),   'prefs');
		
		roiPerimeter = getappdata(hAxis, 'roiPerimeter');
		roiArea      = getappdata(hAxis, 'roiArea');
		roiData      = getappdata(hAxis, 'roiData');
		roiStd       = getappdata(hAxis, 'roiStd');

		% If we need to create new plots, need to ensure that the axis is the current one
		hFig         = get(hAxis, 'Parent');
		oldAxis      = get(hFig,  'CurrentAxes');

		% Escape if invalid iRoi
		if currentTime > numel(roi{currentSlice})
			return
		end
		
		if iRoi > numel(roi{currentSlice}{currentTime})
			return
		end

		% Create plot objects if necessary
		if verLessThan('matlab','8.1.0')
			if isempty(roiHandles)
				roiHandles       = struct('line', [],  'vertices', [],  'center', [],  'label', []);
			end
			if iRoi > numel(roiHandles)
				roiHandles(iRoi) = struct('line', [],  'vertices', [],  'center', [],  'label', []);
			end
		else
			% gobjects are mandatory in R2014b+, but not available until R2013a (8.0.0)
			if isempty(roiHandles)
				roiHandles       = struct('line', gobjects(0,1),  'vertices', gobjects(0,1),  'center', gobjects(0,1),  'label', gobjects(0,1));
			end
			if iRoi > numel(roiHandles)
				roiHandles(iRoi) = struct('line', gobjects(0,1),  'vertices', gobjects(0,1),  'center', gobjects(0,1),  'label', gobjects(0,1));
			end
		end

		% If it's a deleted ROI, hide the plot handles with dummy values
		if isempty(roi{currentSlice}{currentTime}{iRoi})
			if ~isempty(roiHandles(iRoi).line)
				set(roiHandles(iRoi).line,     'XData', [],   'YData', []);
			end
			if ~isempty(roiHandles(iRoi).vertices)
				set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
			end
			if ~isempty(roiHandles(iRoi).center)
				set(roiHandles(iRoi).center,   'XData', [],   'YData', []);
			end
			if ~isempty(roiHandles(iRoi).label)
				set(roiHandles(iRoi).label,    'Position', [0 0 0]);
			end
			return
		end

		% Hide vertices, for the cases where the new ROI has less vertices than the last ROI displayed
		% Bit of a workaround until we can properly hide vertices
% 		if ~isempty(roiHandles(iRoi).vertices) && strcmp(get(roiHandles(iRoi).vertices(1), 'Visible'), 'on')
% 			set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
% 		end
		
		if ~isempty(roiHandles(iRoi).vertices)
			if (verLessThan('matlab','8.4.0') && strcmp(get(roiHandles(iRoi).vertices(1), 'Visible'), 'on')) || (~verLessThan('matlab','8.4.0') && strcmp(roiHandles(iRoi).vertices(1).Visible, 'on'))
				set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
			end
		end
% 		for i = 1:numel(roiHandles(iRoi).vertices)
% 			roiHandles(iRoi).vertices(i).XData = [];
% 			roiHandles(iRoi).vertices(i).YData = [];
% 		end

		% Create plot objects if necessary
		% XLimInclude/YLimInclude are undocumented and used to speed up performance: http://undocumentedmatlab.com/blog/plot-liminclude-properties/
		roiColor = prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:);
		if isempty(roiHandles(iRoi).line)
% 			if (oldAxis ~= hAxis), set(hFig, 'CurrentAxes', hAxis); end
			roiHandles(iRoi).line   = plot(hAxis, 0, 0, ...
		                                 'Color',         roiColor, ...
		                                 'LineStyle',     '-', ...
		                                 'Marker',        'none', ...
		                                 'XLimInclude',   'off', ...
		                                 'YLimInclude',   'off', ...
																		 'ButtonDownFcn', @AddVertexFromEdge);
		end

		% We now do this early on
% 		% Sanity check to avoid plotting way too many points
% 		if any(cellfun(@(x) numel(x), roi{currentSlice}{currentTime}) > 200)
% 			prefs.disableVertices = true;
% 		end
		
		if ~prefs.disableVertices
			for iVertex = 1:size(roi{currentSlice}{currentTime}{iRoi},1)
				if iVertex > length(roiHandles(iRoi).vertices)
	% 				if (oldAxis ~= hAxis), set(hFig, 'CurrentAxes', hAxis); end
					roiHandles(iRoi).vertices(iVertex) = plot(hAxis, 0, 0, ...
					                                          'Color',         roiColor, ...
					                                          'LineStyle',     'none', ...
					                                          'Marker',        '.', ...
					                                          'XLimInclude',   'off', ...
					                                          'YLimInclude',   'off', ...
					                                          'ButtonDownFcn', @VertexButtonDown);
	% 				                                          'MarkerSize',    6, ...
	
					if ismac
						roiHandles(iRoi).vertices(iVertex).MarkerSize = 12;
					end
				end
			end
		end

		if isempty(roiHandles(iRoi).center)
% 			if (oldAxis ~= hAxis), set(hFig, 'CurrentAxes', hAxis); end
			roiHandles(iRoi).center = plot(hAxis, 0, 0, ...
		                                 'Color',      roiColor, ...
		                                 'LineStyle',  'none', ...
		                                 'Marker',     '+', ...
		                                 'XLimInclude',   'off', ...
		                                 'YLimInclude',   'off', ...
		                                 'ButtonDownFcn', @StartRoiMove);
% 		                                 'MarkerSize', 6, ...
		end

		if ~all(prefs.PixelSpacing == 1)
			strUnits = 'mm';
		else
			strUnits = 'px';
		end

% 		if ~prefs.hideRoiLabels && prefs.showRoiInfo
		if any(strcmp(prefs.roiLabelMode, {'area', 'both'}))
			% Perimeter and area measurements
			if (iRoi <= numel(roiPerimeter{currentSlice}{currentTime})) && ~isnan(roiPerimeter{currentSlice}{currentTime}{iRoi})
				strPerimeter = sprintf('%1.1f %s, ', roiPerimeter{currentSlice}{currentTime}{iRoi}, strUnits);
			else
				strPerimeter = '';
			end

			if (iRoi <= numel(roiPerimeter{currentSlice}{currentTime})) && ~isnan(roiArea{currentSlice}{currentTime}{iRoi}) && (roiArea{currentSlice}{currentTime}{iRoi} ~= 0)
				strArea = sprintf('%1.1f %s^2', roiArea{currentSlice}{currentTime}{iRoi}, strUnits);
			else
				strArea = '';
			end
			if ~isempty(strPerimeter) && isempty(strArea)
% 				strPerimeterArea = cat(2, ' (', strPerimeter(1:end-2), ')');
				strPerimeterArea = strPerimeter(1:end-2);
			elseif isempty(strPerimeter) && isempty(strArea)
				strPerimeterArea = '';
			else
				strPerimeterArea = cat(2, strPerimeter, strArea);
			end
		else
			strPerimeterArea = '';
		end

		if any(strcmp(prefs.roiLabelMode, {'signal', 'both'}))
			% Mean and standard deviations
			if (iRoi <= numel(roiData{currentSlice}{currentTime})) && ~isempty(roiData{currentSlice}{currentTime}{iRoi}) && ~isnan(roiData{currentSlice}{currentTime}{iRoi}) && (roiData{currentSlice}{currentTime}{iRoi} ~= 0)
				if roiStd{currentSlice}{currentTime}{iRoi}
					strSignal = sprintf('%1.1f±%1.1f a.u.', roiData{currentSlice}{currentTime}{iRoi}, roiStd{currentSlice}{currentTime}{iRoi});
				else
					strSignal = sprintf('%1.1f a.u.', roiData{currentSlice}{currentTime}{iRoi});
				end
			else
				strSignal = '';
			end
		else
			strSignal = '';
		end
		
		switch(prefs.roiLabelMode)
			case 'number'
				strLabel = num2str(iRoi);
% 				strLabel = sprintf('%1.1f, %1.1f', roi{currentSlice}{currentTime}{iRoi}); % Debug case to show position of points
			case 'area'
				strLabel = sprintf('%d (%s)', iRoi, strPerimeterArea);
			case 'signal'
				strLabel = sprintf('%d (%s)', iRoi, strSignal);
			case 'both'
				strLabel = sprintf('%d (%s)\n%s', iRoi, strSignal, strPerimeterArea);
		end
% 		if prefs.showRoiInfo
% 			strLabel = sprintf('%d%s', iRoi, strPerimeterArea);
% 		else
% 			strLabel = num2str(iRoi);
% 		end

		if isempty(roiHandles(iRoi).label) && ~strcmp(prefs.roiLabelMode, 'none')
			if (oldAxis ~= hAxis), set(hFig, 'CurrentAxes', hAxis); end
			% TODO/FIXME: why was HitTest set to 'on'?
			roiHandles(iRoi).label  = text(0, 0, strLabel, ...
			                               'Parent',              hAxis, ...
			                               'BackgroundColor',     0.25*[1 1 1], ...
			                               'Color',               roiColor, ...
			                               'HitTest',             'off', ...
			                               'HorizontalAlignment', 'right', ...
			                               'VerticalAlignment',   'bottom', ...
			                               'FontName',            prefs.FontName, ...
			                               'FontSize',            prefs.FontSize, ...
																		 'ButtonDownFcn',       @LabelButtonDown);
			if ~verLessThan('matlab','8.4.0')
				roiHandles(iRoi).label.Margin        = 2;
				roiHandles(iRoi).label.FontSmoothing = 'off';
			end
		elseif ~isempty(roiHandles(iRoi).label) && ~strcmp(prefs.roiLabelMode, 'none')
			set(roiHandles(iRoi).label,  'String', strLabel);
		elseif ~isempty(roiHandles(iRoi).label) && strcmp(prefs.roiLabelMode, 'none')
			set(roiHandles(iRoi).label,  'Position', nan(1,2));
		end

		% Update plot object positions
		if ~strcmp(prefs.interpAlgorithm, 'linear') && any(strcmp(prefs.interpRois, {'display', 'always'}))
			tmpInterpRoi = CubicInterp(roi{currentSlice}{currentTime}{iRoi}, [], prefs.interpAlgorithm);
			set(roiHandles(iRoi).line, 'XData', tmpInterpRoi(:,1), 'YData', tmpInterpRoi(:,2));
		else
			% Line segments must be plotted differently as of R2014b
			if any(isnan(roi{currentSlice}{currentTime}{iRoi}))
				set(roiHandles(iRoi).line, 'XData', roi{currentSlice}{currentTime}{iRoi}(1:end-1,1), ...
				                           'YData', roi{currentSlice}{currentTime}{iRoi}(1:end-1,2));
			else
				set(roiHandles(iRoi).line, 'XData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],1), ...
				                           'YData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],2));
			end
		end

		set(roiHandles(iRoi).center, 'XData',    mean(roi{currentSlice}{currentTime}{iRoi}(:,1)), ...
		                             'YData',    mean(roi{currentSlice}{currentTime}{iRoi}(:,2)));
		if ~strcmp(prefs.roiLabelMode, 'none')
			set(roiHandles(iRoi).label,  'Position', [min(roi{currentSlice}{currentTime}{iRoi}(:,1))-0, ...
			                                          min(roi{currentSlice}{currentTime}{iRoi}(:,2))-0 0]);
		end
		% Don't update if they're not visible anyways
		% FIXME: This will break if someone toggles visible back on
		if ~all(strcmp(get(roiHandles(iRoi).vertices, 'Visible'), 'off')) && ~isempty(roiHandles(iRoi).vertices)
			for iVertex = 1:size(roi{currentSlice}{currentTime}{iRoi},1)
% 				if verLessThan('matlab','8.4.0')
% 					roiHandles(iRoi).vertices(iVertex).XData = roi{currentSlice}{currentTime}{iRoi}(iVertex,1);
% 					roiHandles(iRoi).vertices(iVertex).YData = roi{currentSlice}{currentTime}{iRoi}(iVertex,2);
% 				else
					set(roiHandles(iRoi).vertices(iVertex), 'XData', roi{currentSlice}{currentTime}{iRoi}(iVertex,1), 'YData', roi{currentSlice}{currentTime}{iRoi}(iVertex,2));
% 				end
			end
		end
		% This should be faster, but it breaks when finishing an ROI for some reason??
		% Maybe because it sets the XData/YData too quickly, and thus the right
		% click also activates its click-down handler
% 		set(roiHandles(iRoi).vertices, 'XData', roi{currentSlice}{currentTime}{iRoi}(:,1), 'YData', roi{currentSlice}{currentTime}{iRoi}(:,2))

		% Hide vertices, for the cases where the new ROI has less vertices than the last ROI displayed
		for iVertex = size(roi{currentSlice}{currentTime}{iRoi},1)+1:numel(roiHandles(iRoi).vertices)
			set(roiHandles(iRoi).vertices(iVertex), 'XData', [],   'YData', []);
		end

		% Do black and white ROIs for frozen colormap subplot (as they are likely colour pixel maps)
		imgData = getappdata(hAxis, 'imgData');
		
		% prefs.blackWhiteRoi controls ROI colors
		bDoBwRois = false;
		if ~islogical(prefs.blackWhiteRoi)
			indCurrAxes = find(getappdata(hFig, 'hAxes') == hAxis);
			if any(prefs.blackWhiteRoi == indCurrAxes)
				bDoBwRois = true;
			end
		elseif prefs.blackWhiteRoi
			bDoBwRois = true;
		end

		if (ndims(imgData{currentSlice}{currentTime}) == 3) || bDoBwRois
			set(roiHandles(iRoi).line,     'Color',           'k')
			set(roiHandles(iRoi).vertices, 'MarkerEdgeColor', 'w')
		end

		setappdata(hAxis, 'roiHandles', roiHandles);

% 		% Reset CurrentAxes if modified
% 		if oldAxis ~= get(hFig,  'CurrentAxes')
% 			set(hFig, 'CurrentAxes', oldAxis);
% 		end
	end

  % ------------------------------------------------------------
	function UpdateRoiData(src, iRoi, aTimes)
	% - Update mask, data, and plots for a given ROI
	% - Note: UpdateRoiData will call UpdateRoiHandles, but not vice versa!
	% - Updates for CurrentAxes if src is hFig, and hAxis if src is an axis

		% Allow for different arguments
		if strcmp(get(src, 'Type'), 'axes')
			hFig  = get(src, 'Parent');
			hAxis = src;
		elseif strcmp(get(src, 'Type'), 'figure')
			hFig  = src;
			hAxis = get(src, 'CurrentAxes');
		end

		UpdateRoiHandles(hAxis, iRoi);

		% If currently drawing an ROI, don't update data on the fly.  In theory, we
		% might want to allow for updating *other* ROIs, but there's no means current
		% reason to update another ROI while drawing
		if ~isempty(getappdata(hAxis, 'inProgressRoi'))
			return
		end

		prefs        = getappdata(hFig,  'prefs');
		roi          = getappdata(hAxis, 'roi');
		roiMask      = getappdata(hAxis, 'roiMask');
		roiData      = getappdata(hAxis, 'roiData');
		roiStd       = getappdata(hAxis, 'roiStd');
		roiPerimeter = getappdata(hAxis, 'roiPerimeter');
		roiArea      = getappdata(hAxis, 'roiArea');
		roiProfile   = getappdata(hAxis, 'roiProfile');

		imgData      = getappdata(hAxis, 'imgData');
		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');

		imgDataIndexed = getappdata(hAxis, 'imgDataIndexed');
		if ~isempty(imgDataIndexed)
			imgData = imgDataIndexed;
		end

		% Update all slices at once (needed if plotting data for all slices at once,
		% but slow if enabled)
		if prefs.plotShowAllSlices
			aUpdateSlices = 1:numel(imgData);
		else
			aUpdateSlices = currentSlice;
		end

for currentSlice = aUpdateSlices

		if ~exist('aTimes', 'var')
% 			aTimes = 1:numel(roi{currentSlice});
			aTimes = 1:numel(imgData{currentSlice});
		end
		
		badInds = (aTimes > numel(imgData{currentSlice})) | (aTimes < 1);
		if any(badInds)
			aTimes(badInds) = [];
		end

		for iTime = aTimes
			% If staticRoi is off, not all times may have this ROI
			if (iRoi > numel(roi{currentSlice}{iTime})) || isempty(roi{currentSlice}{iTime}{iRoi}) || isempty(imgData{currentSlice}{iTime})
				continue
			end
			if size(roi{currentSlice}{iTime}{iRoi},1) == 1
				% Single point ROIs have a special case
				roiMask{currentSlice}{iTime}{iRoi} = false(size(imgData{currentSlice}{iTime},1), size(imgData{currentSlice}{iTime},2));
% 				roiMask{currentSlice}{iTime}{iRoi}(floor(roi{currentSlice}{iTime}{iRoi}(1,2)), floor(roi{currentSlice}{iTime}{iRoi}(1,1))) = true;
				
				% Round the coordinate, but clip to the sides
				x = min([max([round(roi{currentSlice}{iTime}{iRoi}(1,2)) 0]) size(imgData{currentSlice}{iTime},1)]);
				y = min([max([round(roi{currentSlice}{iTime}{iRoi}(1,1)) 0]) size(imgData{currentSlice}{iTime},2)]);
				roiMask{currentSlice}{iTime}{iRoi}(x,y) = true;
			elseif any(isnan(roi{currentSlice}{iTime}{iRoi}))
				% Line: No data
				roiMask{currentSlice}{iTime}{iRoi} = false(size(imgData{currentSlice}{iTime},1), size(imgData{currentSlice}{iTime},2));
			else
				roiMask{currentSlice}{iTime}{iRoi} = poly2mask(roi{currentSlice}{iTime}{iRoi}(:,1), roi{currentSlice}{iTime}{iRoi}(:,2), size(imgData{currentSlice}{iTime},1), size(imgData{currentSlice}{iTime},2));
			end

			if exist('nanmean.m', 'file')
				roiData{currentSlice}{iTime}{iRoi} = nanmean(double(imgData{currentSlice}{iTime}(roiMask{currentSlice}{iTime}{iRoi})));
			else
				roiData{currentSlice}{iTime}{iRoi} = mean(double(imgData{currentSlice}{iTime}(roiMask{currentSlice}{iTime}{iRoi})));
			end
			
			if exist('nanstd.m', 'file')
				roiStd{ currentSlice}{iTime}{iRoi} = nanstd( double(imgData{currentSlice}{iTime}(roiMask{currentSlice}{iTime}{iRoi})));
			else
				roiStd{ currentSlice}{iTime}{iRoi} = std( double(imgData{currentSlice}{iTime}(roiMask{currentSlice}{iTime}{iRoi})));
			end

			% --- Profile tool for lines ---
			if any(isnan(roi{currentSlice}{iTime}{iRoi}))
				% Hack for RGB data
				tImg = imgData{currentSlice}{iTime};
				if ~ismatrix(tImg)
					tImg = mean(tImg,3);
				end
				
				% Hack for multi-segment lines in a single "ROI"
				if any(isnan(roi{currentSlice}{iTime}{iRoi}(1:end-1,:)))
					roiProfile{currentSlice}{iTime}{iRoi} = nan;
				else
					roiProfile{currentSlice}{iTime}{iRoi} = improfile(tImg, roi{currentSlice}{iTime}{iRoi}(1:end-1,1), roi{currentSlice}{iTime}{iRoi}(1:end-1,2), 'nearest');
				end
			end

			% --- Calculate perimeter ---
			tmpRoi = roi{currentSlice}{iTime}{iRoi};

			if isfield(prefs, 'PixelSpacing') && ~isempty(prefs.PixelSpacing)
				tmpRoi(:,1) = tmpRoi(:,1) .* prefs.PixelSpacing(2);
				tmpRoi(:,2) = tmpRoi(:,2) .* prefs.PixelSpacing(1);
			end

			% Close/unclose the polygon if necessary
			if (sum(abs(diff(tmpRoi([1 end],:)))) ~= 0) && (size(tmpRoi,1) > 2)
				tmpRoi = tmpRoi([1:end 1],:);
			elseif (sum(abs(diff(tmpRoi([1 end],:)))) < 1e-9) && (size(tmpRoi,1) > 1) && (size(tmpRoi,1) < 4)
				tmpRoi = tmpRoi(1:end-1,:);
			end

			if size(tmpRoi,1) == 1
				roiPerimeter{currentSlice}{iTime}{iRoi} = nan;
			else
% 				roiPerimeter{currentSlice}{iTime}{iRoi} = sum(sqrt(sum(diff(tmpRoi([1:end 1],:)).^2,2)));
				% Special case for line segments
				if all(isnan(roi{currentSlice}{iTime}{iRoi}(end,:)))
					roiPerimeter{currentSlice}{iTime}{iRoi} = sum(sqrt(sum(diff(tmpRoi(1:end-2,:)).^2,2)));
				else
					roiPerimeter{currentSlice}{iTime}{iRoi} = sum(sqrt(sum(diff(tmpRoi).^2,2)));
				end
			end

			roiArea{currentSlice}{iTime}{iRoi} = polyarea(tmpRoi(:,1), tmpRoi(:,2));
		end
end
		setappdata(hAxis, 'roiMask',      roiMask);
		setappdata(hAxis, 'roiData',      roiData);
		setappdata(hAxis, 'roiStd',       roiStd);
		setappdata(hAxis, 'roiPerimeter', roiPerimeter);
		setappdata(hAxis, 'roiArea',      roiArea);
		setappdata(hAxis, 'roiProfile',   roiProfile);

		UpdateRoiFigure(hFig)
	end

	%% ------------------------------------------------------------
	function LabelButtonDown(src, evt)
		% Code for these conditions is from VertexButtonDown()
		if strcmp(get(src, 'Type'), 'figure')
			hFig = src;
			hObj = get(hFig, 'CurrentObject');
		else
			hFig = get(get(src, 'Parent'), 'Parent');
			hObj = src;
		end

		% Don't allow this when nudging, as it gets in the way and we already have a way of adjusting the ROI
		if any(strcmp(getappdata(hFig, 'strMode'), {'Nudge', 'Circle'}))
			return
		end

		% Also not allowed while drawing
		if ~isempty(getappdata(get(hFig, 'CurrentAxes'), 'inProgressRoi'))
			return
		end

		if strcmp(get(hFig, 'SelectionType'), 'extend')
			% Shift-click to toggle background color
			if strcmp(get(hObj, 'BackgroundColor'), 'none')
				set(hObj, 'BackgroundColor', [0 0 0])
			else
				set(hObj, 'BackgroundColor', 'none')
			end
		else
			% Store the delta between the clicked point and the stored label position
			hAxis = get(hFig, 'CurrentAxes');
			LabelPosition = get(hObj,  'Position');
			CurrentPoint  = get(hAxis, 'CurrentPoint');
			setappdata(hFig, 'LabelDelta', LabelPosition(1,1:2) - CurrentPoint(1,1:2));

			set(hFig, 'WindowButtonMotionFcn', @LabelMove);
			set(hFig, 'WindowButtonUpFcn',     @EndLabelMove);
		end
	end

  % ------------------------------------------------------------
	function EndLabelMove(src, evt)
		hFig = src;

		set(hFig, 'WindowButtonMotionFcn', '');
		set(hFig, 'WindowButtonUpFcn',     '');
	end

  % ------------------------------------------------------------
	function LabelMove(hFig, evt)
		hAxis        = get(hFig,        'CurrentAxes');
		hObj         = get(hFig,        'CurrentObject');
		LabelDelta   = getappdata(hFig, 'LabelDelta');
		CurrentPoint = get(hAxis,       'CurrentPoint');
		
		% LabelDelta accounts for difference between clicked point and the stored label position
		set(hObj, 'Position',  [CurrentPoint(1,1:2)+LabelDelta 0])
	end

	%% ------------------------------------------------------------
	function CreateRoiFigure(hMainFig)
% 		setappdata(hFig, 'roiPlot',     struct('hFig', [], 'hAxes', [])
		roiPlot   = getappdata(hMainFig, 'roiPlot');
		mainPos   = get(       hMainFig, 'Position');  % [left, bottom, width, height]
		hMainAxes = getappdata(hMainFig, 'hAxes');
		prefs     = getappdata(hMainFig, 'prefs');

		if ~isempty(roiPlot.hFig) && ishandle(roiPlot.hFig)
			return
		end

		newFigPos = mainPos + [mainPos(3) 0 0 0];  % Right of main figure

		% Check to see if the new figure will be pushed off the edge of the screen
		% If so, push the main figure left as much as possible
		previousUnits = get(0, 'Units');
		set(0, 'Units', 'pixels');
		screenPos = get(0, 'ScreenSize');  % [left bottom width height]
		set(0, 'Units', previousUnits);

		if (newFigPos(1) + newFigPos(3)) > screenPos(3)
			mainPos(1) = max([0 screenPos(3) - 2*mainPos(3)]);
			newFigPos(1) = mainPos(1) + mainPos(3);
		end

		% Create it invisible first to minimize the focus stealing time
		hFig = figure;
		set(hFig, 'Toolbar',         'none');
		set(hFig, 'Color',           0.75*[1 1 1]);
		set(hFig, 'NumberTitle',     'off');
		set(hFig, 'KeyPressFcn',     @CatchKeyPressRoiPlot);
		set(hFig, 'CloseRequestFcn', @CloseRoiFigure);
		setappdata(hFig, 'hMainFig', hMainFig);

		set(hMainFig, 'Position', mainPos);
		set(hFig,     'Position', newFigPos);

		% Create UI tabs
		if verLessThan('matlab','8.4.0')
			ws = warning('off', 'MATLAB:uitabgroup:OldVersion');
		end

		hTabGroup = uitabgroup('Parent', hFig);
		setappdata(hFig, 'hTabGroup', hTabGroup);

		if ~verLessThan('matlab','8.4.0')
			hTabGroup.SelectionChangedFcn = @SwitchRoiFigTab;
		else
			set(hTabGroup, 'SelectionChangeFcn', @SwitchRoiFigTab);
			warning(ws);
		end
		
		htTableMean = uitab('Parent', hTabGroup, 'Title', 'Mean');     setappdata(hMainFig, 'htTableMean', htTableMean);
		htTableStd  = uitab('Parent', hTabGroup, 'Title', 'Std Dev');  setappdata(hMainFig, 'htTableStd',  htTableStd);
		htPlot      = uitab('Parent', hTabGroup, 'Title', 'Plot');     setappdata(hFig,     'htPlot',      htPlot);

		if prefs.showMMode
			htMMode   = uitab('Parent', hTabGroup, 'Title', 'M-Mode');    setappdata(hFig, 'htMMode',   htMMode);
		end

		% Default to the plot
		set(hTabGroup,'SelectedTab', htPlot);

		% Remember, hAxes is the transpose of the subplot layout
		nCols = size(hMainAxes,1);
		nRows = size(hMainAxes,2);

		hAxesPlot   = zeros(size(hMainAxes));
		hAxesMMode  = zeros(size(hMainAxes));
		hTableMean  = zeros(size(hMainAxes));
		hTableStd   = zeros(size(hMainAxes));
		hImageMMode = zeros(size(hMainAxes));
		for iPlot = 1:numel(hMainAxes)
			if hMainAxes(iPlot) == 0
				continue
			end

			% Create raw data plot axes
			hAxesPlot(iPlot) = subplottight(nRows, nCols, iPlot, 0.85, 'Parent', htPlot);
			grid on
			if verLessThan('matlab','8.4.0')
				set(hAxesPlot(iPlot), 'GridLineStyle', '--')
			end
			set(hAxesPlot(iPlot), 'Color',  0.7*[1 1 1])%'none')
			set(hAxesPlot(iPlot), 'XColor', 'k') % 'w')
			set(hAxesPlot(iPlot), 'YColor', 'k') % 'w')
			axis off
			xlabel('Time'), ylabel('Signal Intensity')

			if prefs.showMMode
				% Create m-mode axes
				hAxesMMode(iPlot) = subplottight(nRows, nCols, iPlot, 0.95, 'Parent', htMMode);
				hImageMMode(iPlot) = imagesc([]);
				axis tight
				set(hAxesMMode(iPlot), 'XTick', [])
				set(hAxesMMode(iPlot), 'YTick', [])
				set(hAxesMMode(iPlot), 'XLimMode', 'auto')
				set(hAxesMMode(iPlot), 'YLimMode', 'auto')
% 				xlabel('Time --->'), ylabel('<--- Spatial Profile')
				xlabel('Time'), ylabel('Spatial Profile')
% 				colormap(gray)
				colormap(get(hMainFig, 'colormap'))
			end
			
			% Raw data table
			% Need to calculate the position (equivalent to subplottight's calculation)
			width  = 1 / nCols;
			height = 1 / nRows;

			fillRatio = 0.98;
			gapX = width  * (1-fillRatio);
			gapY = height * (1-fillRatio);
			[iCol, iRow] = ind2sub([nCols nRows], iPlot);

			left   = (iCol - 1)     * width;
			bottom = (nRows - iRow) * height;

			position = [left + gapX, bottom + gapY, width - 2*gapX, height - 2*gapY];

			hTableMean(iPlot) = uitable(htTableMean, ...
            'Data',         [], ...
            'ColumnName',   {'test'}, ...
            'RowName',      {'test row'}, ...
            'Units',        'normalized', ...
            'Position',     position, ...
            'ColumnWidth',  'auto', ...
            'ColumnFormat', {'shortg'});
%             'ColumnWidth',  {100});

			hTableStd(iPlot) = uitable(htTableStd, ...
            'Data',         [], ...
            'ColumnName',   {'test'}, ...
            'RowName',      {'test row'}, ...
            'Units',        'normalized', ...
            'Position',     position, ...
            'ColumnWidth',  'auto', ...
            'ColumnFormat', {'shortg'});
%             'ColumnWidth',  {100});
		end

		roiPlot.hFig  = hFig;
		roiPlot.hAxes = hAxesPlot;
		setappdata(hMainFig, 'roiPlot', roiPlot);

		% This seems like a poor way of doing it, considering roiPlot and the fact
		% that we'll need htRawData later too...
		setappdata(hMainFig, 'hTabGroup',   hTabGroup);
		setappdata(hMainFig, 'hTableMean',  hTableMean);
		setappdata(hMainFig, 'hTableStd',   hTableStd);
		if prefs.showMMode
			setappdata(hMainFig, 'hAxesMMode',  hAxesMMode);
			setappdata(hMainFig, 'hImageMMode', hImageMMode);
		end
		drawnow % Force the new figure to be created so that we can restore focus to the main figure
% 		set(hFig, 'Visible', 'on')
		figure(hMainFig)
	end

	%% ---------------------------------------------------------------------------
	function SwitchRoiFigTab(src, evt)
	% Update the data for the tab that we switch to
		if verLessThan('matlab','8.4.0')
			hTabGroupChildren = get(src, 'Children');
			strSelectedTabTitle = get(hTabGroupChildren(evt.NewValue), 'Title');
			src = handle(src);
		else
			strSelectedTabTitle = evt.NewValue.Title;
		end

		switch strSelectedTabTitle
			case 'M-Mode'
				hMainFig = getappdata(src.Parent, 'hMainFig');
				fcnHandles = getappdata(hMainFig, 'fcnHandles');
				fcnHandles.UpdateMMode(hMainFig);
			case 'Mean'
				hMainFig = getappdata(src.Parent, 'hMainFig');
				fcnHandles = getappdata(hMainFig, 'fcnHandles');
				fcnHandles.UpdateRawTable(hMainFig, strSelectedTabTitle);
			case 'Std Dev'
				hMainFig = getappdata(src.Parent, 'hMainFig');
				fcnHandles = getappdata(hMainFig, 'fcnHandles');
				fcnHandles.UpdateRawTable(hMainFig, strSelectedTabTitle);
			otherwise
				return
		end
end
	
	%% ------------------------------------------------------------
	function CloseRoiFigure(src, evt)
		hMainFig = getappdata(src, 'hMainFig');
		if ishandle(hMainFig)
			setappdata(hMainFig, 'roiPlot', struct('hFig', [], 'hAxes', []));

			% Disable the ROI figure for the future too
			prefs = getappdata(hMainFig, 'prefs');
			prefs.suppressRoiPlot = true;
			setappdata(hMainFig, 'prefs', prefs);
		end
		delete(src);
	end

	function ToggleRoiFigure(hFig)
		% Toggles whether the ROI plot/graph is displayed or not
		prefs = getappdata(hFig, 'prefs');

		if prefs.suppressRoiPlot
			prefs.suppressRoiPlot = false;
			setappdata(hFig, 'prefs', prefs);
			CreateRoiFigure(hFig);
			UpdateRoiFigure(hFig);
		else
			roiPlot = getappdata(hFig, 'roiPlot');
			if ~isempty(roiPlot) && isfield(roiPlot, 'hFig') && ~isempty(roiPlot.hFig) && ishandle(roiPlot.hFig)
				CloseRoiFigure(roiPlot.hFig); % This also toggles prefs.suppressRoiPlot
				return
			end
		end
	end
	
	%% ------------------------------------------------------------
	function UpdateRoiFigure(hMainFig)
		prefs = getappdata(hMainFig, 'prefs');
		hTabGroup = getappdata(hMainFig, 'hTabGroup');
		if prefs.suppressRoiPlot
			return
		end

		if prefs.plotShowAllSlices
			UpdateRoiFigureAllSlices(hMainFig);
			return
		end

		if prefs.showMMode && ~isempty(hTabGroup) && ishandle(hTabGroup)
			if verLessThan('matlab','8.4.0')
				hTabGroupChildren = get(hTabGroup, 'Children');
				strSelectedTabTitle = get(hTabGroupChildren(get(hTabGroup, 'SelectedIndex')), 'Title');
			else
				strSelectedTabTitle = hTabGroup.SelectedTab.Title;
			end
			if strcmp(strSelectedTabTitle, 'M-Mode')
				UpdateMMode(hMainFig);
			end
		end

		if ~isempty(hTabGroup) && ishandle(hTabGroup)
			if verLessThan('matlab','8.4.0')
				hTabGroupChildren = get(hTabGroup, 'Children');
				strSelectedTabTitle = get(hTabGroupChildren(get(hTabGroup, 'SelectedIndex')), 'Title');
			else
				strSelectedTabTitle = hTabGroup.SelectedTab.Title;
			end
			if strcmp(strSelectedTabTitle, 'Mean')
				UpdateRawTable(hMainFig, strSelectedTabTitle);
			end
		end

		roiPlot      = getappdata(hMainFig, 'roiPlot');
		hMainAxes    = getappdata(hMainFig, 'hAxes');
		prefs        = getappdata(hMainFig, 'prefs');

		% Too slow... if ROI figure is still active from last time updated, this fails
		hAxis        = get(hMainFig, 'CurrentAxes');

		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');

		% If no valid ROI plot, create if necessary, otherwise quit
		if isempty(roiPlot.hFig) || ~ishandle(roiPlot.hFig)
			bNoRoiData = true;
			for jAxis = reshape(hMainAxes(hMainAxes ~= 0), 1, [])
				roiData = getappdata(jAxis, 'roiData');
				if (numel(roiData{currentSlice}) >= currentTime) && ~isempty(roiData{currentSlice}{currentTime})
					bNoRoiData = false;
				end
			end
			if bNoRoiData
				return
			end
			CreateRoiFigure(hMainFig);
			roiPlot = getappdata(hMainFig, 'roiPlot');
		end

		for iRow = 1:size(roiPlot.hAxes,1)
			for iCol = 1:size(roiPlot.hAxes,2)
				if hMainAxes(iRow,iCol) == 0
					continue
				end

				hPlots       = getappdata(roiPlot.hAxes(iRow,iCol), 'hPlots');
				hPlotStd     = getappdata(roiPlot.hAxes(iRow,iCol), 'hPlotStd');
				roi          = getappdata(hMainAxes(iRow,iCol),     'roi');
				roiData      = getappdata(hMainAxes(iRow,iCol),     'roiData');
				roiStd       = getappdata(hMainAxes(iRow,iCol),     'roiStd');
				roiProfile   = getappdata(hMainAxes(iRow,iCol),     'roiProfile');
				roiPerimeter = getappdata(hMainAxes(iRow,iCol),     'roiPerimeter');

				if isempty(roiStd)
					roiStd = cell(1,0);
				end

				% Bypass this if there's no ROI data for this current time (happens if
				% subplots have differet number of frames)
				if (numel(roiData{currentSlice}) < currentTime) || (currentTime == 0)
					continue
				end

				% Remove all plots if no ROIs
				if isempty(roiData{currentSlice}{currentTime})
					if ~isempty(hPlots)
						delete(hPlots(~isnan(hPlots)))
						set(roiPlot.hAxes(iRow,iCol), 'Visible', 'off')
						setappdata(roiPlot.hAxes(iRow,iCol), 'hPlots', []);
					end

					if prefs.showStdPlot
						if ~isempty(hPlotStd)
							if ~iscell(hPlotStd)
								disp('hPlotStd isn''t a cell')
							else
% 								tmp = cell2mat(hPlotStd);
% 								delete(tmp(~isnan(tmp)))
								for i = 1:numel(hPlotStd)
									delete(hPlotStd{i})
								end
								setappdata(roiPlot.hAxes(iRow,iCol), 'hPlotStd', cell(1,0));
							end
						end
					end
					continue
				end

				% Turn on visibility if it's not already so
				if strcmp(get(roiPlot.hAxes(iRow,iCol), 'Visible'), 'off')
					set(roiPlot.hAxes(iRow,iCol), 'Visible', 'on')
				end

				% If the plots don't match the data, redo them all
				if numel(roiData{currentSlice}{currentTime}) ~= numel(hPlots)
					delete(hPlots  ),           hPlots   = [];
					for i = 1:numel(hPlotStd)
						delete(hPlotStd{i})
					end
					hPlotStd = cell(1,0);
					figure(roiPlot.hFig)
					subplot(roiPlot.hAxes(iRow,iCol))
					hold on

					for iRoi = 1:length(roiData{currentSlice}{currentTime})
						hPlots(iRoi) = plot(roiPlot.hAxes(iRow,iCol), 1,1, '.-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
					end

					if prefs.showStdPlot
						for iRoi = 1:length(roiData{currentSlice}{currentTime})
							for iTime = 1:length(roiData{currentSlice})
								hPlotStd{iRoi}(iTime) = plot(roiPlot.hAxes(iRow,iCol), [1 1], [1 1], '-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
							end
						end
					end

					% Decide if we're plotting profiles
					bPlotProfile = any(cellfun(@(x) any(cellfun(@(y) any(isnan(y(:))), x)), roi{1})); % Check for a NaN in any roi for all times in current slice

					if bPlotProfile
% 						set(roiPlot.hAxes(iRow,iCol), 'XLim', [0.75 100.25])
						if ~all(prefs.PixelSpacing == 1)
							strUnits = 'mm';
						else
							strUnits = 'px';
						end
% 						maxProfileLength = max(cellfun(@(x) max(cell2mat(x)), roiPerimeter{currentSlice}{currentTime}));
						maxProfileLength = max(cell2mat(roiPerimeter{currentSlice}{currentTime}));
						set(roiPlot.hAxes(iRow,iCol), 'XLim', [0 maxProfileLength])
						set(get(roiPlot.hAxes(iRow,iCol), 'XLabel'), 'String', sprintf('Spatial Profile (%s)', strUnits))
					else
						% Find some non-empty index to get the proper length
						set(roiPlot.hAxes(iRow,iCol), 'XLim', [0.75 length(roiData{currentSlice})+.25])
					end
					set(roiPlot.hAxes(iRow,iCol), 'YLimMode', 'auto')
					figure(hMainFig)
				end

				% Update plots for ROIs
				for iRoi = 1:length(roiData{currentSlice}{currentTime})
					if isempty(roiData{currentSlice}{currentTime}{iRoi})
						if ~isnan(hPlots(iRoi))
							delete(hPlots(iRoi));
							hPlots(iRoi) = nan;
						end

						if prefs.showStdPlot
							if ~isempty(hPlotStd{iRoi})
								delete(hPlotStd{iRoi});
								hPlotStd{iRoi} = [];
							end
						end
						continue
					end

					if isnan(hPlots(iRoi))
						hPlots(iRoi) = plot(roiPlot.hAxes(iRow,iCol), 1,1, '.-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
					end

					if prefs.showStdPlot
						if isempty(hPlotStd) || isempty(hPlotStd{iRoi})
							for iTime = 1:length(roiData{currentSlice})
								hPlotStd{iRoi}(iTime) = plot(roiPlot.hAxes(iRow,iCol), [1 1],[1 1], '-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
							end
						end
					end

					dataToPlot = nan(1, numel(roiData{currentSlice}));
					stdToPlot  = nan(1, numel(roiData{currentSlice}));
					for iTime = 1:numel(roiData{currentSlice})
						if (numel(roiData{currentSlice}{iTime}) >= iRoi) && ~isempty(roiData{currentSlice}{iTime}{iRoi})
							dataToPlot(iTime) = roiData{currentSlice}{iTime}{iRoi};
							stdToPlot( iTime) = roiStd{currentSlice}{iTime}{iRoi};
						end

						if prefs.showStdPlot && (stdToPlot(iTime) ~= 0)
							set(hPlotStd{iRoi}(iTime), 'XData', iTime*[1 1])
							set(hPlotStd{iRoi}(iTime), 'YData', dataToPlot(iTime) + stdToPlot(iTime)*[-1 1])
						else
							if ~isempty(hPlotStd)
								set(hPlotStd{iRoi}(iTime), 'XData', nan(1,2))
								set(hPlotStd{iRoi}(iTime), 'YData', nan(1,2))
							end
						end							
					end

					% Show profile data for lines (with NaNs)
					if (iRoi <= numel(roi{currentSlice}{currentTime})) && any(isnan(roi{currentSlice}{currentTime}{iRoi}(:)))
% 						set(hPlots(iRoi), 'XData', 1:numel(roiProfile{currentSlice}{currentTime}{iRoi}))
% 						set(hPlots(iRoi), 'YData', roiProfile{currentSlice}{currentTime}{iRoi})
% 						set(hPlots(iRoi), 'XData', linspace(1,100, numel(roiProfile{currentSlice}{currentTime}{iRoi})))
						set(hPlots(iRoi), 'XData', linspace(0,roiPerimeter{currentSlice}{currentTime}{iRoi}, numel(roiProfile{currentSlice}{currentTime}{iRoi})))
						set(hPlots(iRoi), 'YData', roiProfile{currentSlice}{currentTime}{iRoi})
						
						% Need to update the XLim here because the profile might change length as it's adjusted
% 						maxProfileLength = max(cellfun(@(x) max(cell2mat(x)), roiPerimeter{currentSlice}));
						maxProfileLength = max(cell2mat(roiPerimeter{currentSlice}{currentTime}));
						set(roiPlot.hAxes(iRow,iCol), 'XLim', [0 maxProfileLength])
					else
						set(hPlots(iRoi), 'XData', 1:numel(dataToPlot))
						set(hPlots(iRoi), 'YData', dataToPlot)
					end
				end

				setappdata(roiPlot.hAxes(iRow,iCol), 'hPlots', hPlots);
				setappdata(roiPlot.hAxes(iRow,iCol), 'hPlotStd', hPlotStd);
			end
		end
% 		figure(hMainFig)
	end

	%% ------------------------------------------------------------
	function UpdateRoiFigureAllSlices(hMainFig)
		prefs = getappdata(hMainFig, 'prefs');
		if prefs.suppressRoiPlot
			return
		end

		roiPlot      = getappdata(hMainFig, 'roiPlot');
		hMainAxes    = getappdata(hMainFig, 'hAxes');
		prefs        = getappdata(hMainFig, 'prefs');

		% Too slow... if ROI figure is still active from last time updated, this fails
		hAxis        = get(hMainFig, 'CurrentAxes');

		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');

		% If no valid ROI plot, create if necessary, otherwise quit
		if isempty(roiPlot.hFig) || ~ishandle(roiPlot.hFig)
			bNoRoiData = true;
			for jAxis = reshape(hMainAxes(hMainAxes ~= 0), 1, [])
				roiData = getappdata(jAxis, 'roiData');
				if (numel(roiData{currentSlice}) >= currentTime) && ~isempty(roiData{currentSlice}{currentTime})
					bNoRoiData = false;
				end
			end
			if bNoRoiData
				return
			end
			CreateRoiFigure(hMainFig);
			roiPlot = getappdata(hMainFig, 'roiPlot');
		end

		for iRow = 1:size(roiPlot.hAxes,1)
			for iCol = 1:size(roiPlot.hAxes,2)
				if hMainAxes(iRow,iCol) == 0
					continue
				end

				hPlots     = getappdata(roiPlot.hAxes(iRow,iCol), 'hPlots');
				hPlotStd   = getappdata(roiPlot.hAxes(iRow,iCol), 'hPlotStd');
				roi        = getappdata(hMainAxes(iRow,iCol),     'roi');
				roiData    = getappdata(hMainAxes(iRow,iCol),     'roiData');
				roiStd     = getappdata(hMainAxes(iRow,iCol),     'roiStd');
				roiProfile = getappdata(hMainAxes(iRow,iCol),     'roiProfile');

				if isempty(roiStd)
					roiStd = cell(1,0);
				end

				% Remove all plots if no ROIs
				if isempty(roiData{currentSlice}{currentTime})
					if ~isempty(hPlots)
						delete(hPlots(~isnan(hPlots)))
						set(roiPlot.hAxes(iRow,iCol), 'Visible', 'off')
						setappdata(roiPlot.hAxes(iRow,iCol), 'hPlots', []);
					end

					if prefs.showStdPlot
						if ~isempty(hPlotStd)
							if ~iscell(hPlotStd)
								disp('hPlotStd isn''t a cell')
							else
% 								tmp = cell2mat(hPlotStd);
% 								delete(tmp(~isnan(tmp)))
								for i = 1:numel(hPlotStd)
									delete(hPlotStd{i})
								end
								setappdata(roiPlot.hAxes(iRow,iCol), 'hPlotStd', cell(1,0));
							end
						end
					end
					continue
				end

				% Turn on visibility if it's not already so
				if strcmp(get(roiPlot.hAxes(iRow,iCol), 'Visible'), 'off')
					set(roiPlot.hAxes(iRow,iCol), 'Visible', 'on')
				end

				% If the plots don't match the data, redo them all
% 				if numel(roiData{currentSlice}{currentTime}) ~= numel(hPlots) % 2014-09-16
				if numel(hPlots) ~= numel(roiData)
					delete(hPlots  ),           hPlots   = [];
% 					delete(cell2mat(hPlotStd)), hPlotStd = cell(1,0);
					for i = 1:numel(hPlotStd)
						delete(hPlotStd{i})
					end
					hPlotStd = cell(1,0);

					figure(roiPlot.hFig)
					subplot(roiPlot.hAxes(iRow,iCol))
					hold on

% 					for iRoi = 1:length(roiData{currentSlice}{currentTime}) % 2014-09-16
					for iRoi = 1:length(roiData) % 2014-09-16
						hPlots(iRoi) = plot(roiPlot.hAxes(iRow,iCol), 1,1, '.-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
					end

					if prefs.showStdPlot
% 						for iRoi = 1:length(roiData{currentSlice}{currentTime}) % 2014-09-16
						for iRoi = 1:length(roiData)
							for iTime = 1:length(roiData{currentSlice})
								hPlotStd{iRoi}(iTime) = plot(roiPlot.hAxes(iRow,iCol), [1 1], [1 1], '-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
							end
						end
					end

					% Decide if we're plotting profiles
					bPlotProfile = any(cellfun(@(x) any(cellfun(@(y) any(isnan(y(:))), x)), roi{1})); % Check for a NaN in any roi for all times in current slice

					if bPlotProfile
						set(roiPlot.hAxes(iRow,iCol), 'XLim', [0.75 100.25])
					else
						% Find some non-empty index to get the proper length
						set(roiPlot.hAxes(iRow,iCol), 'XLim', [0.75 length(roiData{currentSlice})+.25])
					end
					set(roiPlot.hAxes(iRow,iCol), 'YLimMode', 'auto')
					figure(hMainFig)
				end

				% Update plots for ROIs
% 				for iRoi = 1:length(roiData{currentSlice}{currentTime}) % 2014-09-16
				iRoi = 1;
				for currentSlice = 1:length(roiData)
					% 2014-09-16 Need another check to make sure the data has been updated
					if numel(roiData{currentSlice}{currentTime}) < iRoi
						continue
					end
					if isempty(roiData{currentSlice}{currentTime}{iRoi})
						if ~isnan(hPlots(currentSlice))
							delete(hPlots(currentSlice));
							hPlots(currentSlice) = nan;
						end

						if prefs.showStdPlot
							if ~isempty(hPlotStd{currentSlice})
								delete(hPlotStd{currentSlice});
								hPlotStd{currentSlice} = [];
							end
						end
						continue
					end

					if isnan(hPlots(currentSlice))
						hPlots(currentSlice) = plot(roiPlot.hAxes(iRow,iCol), 1,1, '.-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
					end

					if prefs.showStdPlot
						if isempty(hPlotStd) || isempty(hPlotStd{iRoi})
							for iTime = 1:length(roiData{currentSlice})
								hPlotStd{currentSlice}(iTime) = plot(roiPlot.hAxes(iRow,iCol), [1 1],[1 1], '-', 'Color', prefs.roiColors(mod(iRoi-1, size(prefs.roiColors,1))+1,:));
							end
						end
					end

					dataToPlot = nan(1, numel(roiData{currentSlice}));
					stdToPlot  = nan(1, numel(roiData{currentSlice}));
					for iTime = 1:numel(roiData{currentSlice})
						if (numel(roiData{currentSlice}{iTime}) >= iRoi) && ~isempty(roiData{currentSlice}{iTime}{iRoi})
							dataToPlot(iTime) = roiData{currentSlice}{iTime}{iRoi};
							stdToPlot( iTime) = roiStd{currentSlice}{iTime}{iRoi};
						end

						if prefs.showStdPlot && (stdToPlot(iTime) ~= 0)
							set(hPlotStd{currentSlice}(iTime), 'XData', iTime*[1 1])
							set(hPlotStd{currentSlice}(iTime), 'YData', dataToPlot(iTime) + stdToPlot(iTime)*[-1 1])
						else
							if ~isempty(hPlotStd)
								set(hPlotStd{currentSlice}(iTime), 'XData', nan(1,2))
								set(hPlotStd{currentSlice}(iTime), 'YData', nan(1,2))
							end
						end
					end

					% Show profile data for lines (with NaNs)
					if (iRoi <= numel(roi{currentSlice}{currentTime})) && any(isnan(roi{currentSlice}{currentTime}{iRoi}(:)))
% 						set(hPlots(iRoi), 'XData', 1:numel(roiProfile{currentSlice}{currentTime}{iRoi}))
% 						set(hPlots(iRoi), 'YData', roiProfile{currentSlice}{currentTime}{iRoi})
						set(hPlots(currentSlice), 'XData', linspace(1,100, numel(roiProfile{currentSlice}{currentTime}{iRoi})))
						set(hPlots(currentSlice), 'YData', roiProfile{currentSlice}{currentTime}{iRoi})
					else
						set(hPlots(currentSlice), 'XData', 1:numel(dataToPlot))
						set(hPlots(currentSlice), 'YData', dataToPlot)
					end
				end

				setappdata(roiPlot.hAxes(iRow,iCol), 'hPlots', hPlots);
				setappdata(roiPlot.hAxes(iRow,iCol), 'hPlotStd', hPlotStd);
			end
		end
% 		figure(hMainFig)
	end

	%% ---------------------------------------------------------------------------
	function UpdateMMode(hMainFig)
		roiPlot      = getappdata(hMainFig, 'roiPlot');
		hMainAxes    = getappdata(hMainFig, 'hAxes');
		prefs        = getappdata(hMainFig, 'prefs');

		hAxesMMode =  getappdata(hMainFig, 'hAxesMMode');
		hImageMMode = getappdata(hMainFig, 'hImageMMode');

		% Too slow... if ROI figure is still active from last time updated, this fails
		hAxis        = get(hMainFig, 'CurrentAxes');

		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');

		% If no valid ROI plot, create if necessary, otherwise quit
		if isempty(roiPlot.hFig) || ~ishandle(roiPlot.hFig)
			bNoRoiData = true;
			for jAxis = reshape(hMainAxes(hMainAxes ~= 0), 1, [])
				roiData = getappdata(jAxis, 'roiData');
				if (numel(roiData{currentSlice}) >= currentTime) && ~isempty(roiData{currentSlice}{currentTime})
					bNoRoiData = false;
				end
			end
			if bNoRoiData
				return
			end
			CreateRoiFigure(hMainFig);
			roiPlot = getappdata(hMainFig, 'roiPlot');
			hAxesMMode =  getappdata(hMainFig, 'hAxesMMode');
			hImageMMode = getappdata(hMainFig, 'hImageMMode');
		end

		for iRow = 1:size(roiPlot.hAxes,1)
			for iCol = 1:size(roiPlot.hAxes,2)
				if hMainAxes(iRow,iCol) == 0
					continue
				end

				roi          = getappdata(hMainAxes(iRow,iCol),     'roi');
				roiProfile   = getappdata(hMainAxes(iRow,iCol),     'roiProfile');
				roiPerimeter = getappdata(hMainAxes(iRow,iCol),     'roiPerimeter');

				% Remove image data
				if isempty(roi{currentSlice}{currentTime}) && ~isempty(get(hImageMMode(iRow,iCol), 'CData'))
					set(hImageMMode(iRow,iCol), 'CData', []);
					continue
				end

				% Update m-mode data
				iRoi = 1;  % Can only display one ROI's data

				% Abort if no ROI or the ROI isn't a profile
				if isempty(roi{currentSlice}{1}) || ~any(isnan(roi{currentSlice}{1}{iRoi}(:)))
					continue
				end

				% Somewhat complicated by the fact that there might be empty time images
				indsGood = find(cellfun(@(x) ~isempty(x), roiProfile{currentSlice}));
				imgMMode = nan(numel(roiProfile{currentSlice}{indsGood(1)}{iRoi}), numel(roiProfile{currentSlice}));
				imgMMode(:,indsGood) = cell2mat(cellfun(@(x) x{iRoi}, roiProfile{currentSlice}(indsGood), 'UniformOutput', false));

				set(hImageMMode(iRow,iCol), 'CData', imgMMode)

				% This is not the right way of doing this
% 				% Use the true distance, if available
% 				set(hImageMMode(iRow,iCol), 'YData', linspace(0, roiPerimeter{currentSlice}{1}{iRoi}, numel(roiProfile{currentSlice}{1}{iRoi})))
				
				if verLessThan('matlab','8.4.0')
					% Need to manually update axis limits if setting for the first time
					xlim = get(get(hImageMMode(iRow,iCol), 'Parent'), 'XLim');
					ylim = get(get(hImageMMode(iRow,iCol), 'Parent'), 'YLim');
					if (xlim(1) ~= 0.5) || (diff(xlim) ~= size(imgMMode,2))
						set(get(hImageMMode(iRow,iCol), 'Parent'), 'XLim', [0.5 size(imgMMode,2)+0.5]);
					end
					
					if (ylim(1) ~= 0.5) || (diff(ylim) ~= size(imgMMode,1))
						set(get(hImageMMode(iRow,iCol), 'Parent'), 'YLim', [0.5 size(imgMMode,1)+0.5]);
					end
				end
			end
		end
	end
	
	%% ---------------------------------------------------------------------------
	function UpdateRawTable(hFig, strSelectedTabTitle)
		roiPlot      = getappdata(hFig,  'roiPlot');
		hMainAxes    = getappdata(hFig,  'hAxes');
		hTableMean   = getappdata(hFig,  'hTableMean');
		hTableStd    = getappdata(hFig,  'hTableStd');

		hAxis        = get(hFig,         'CurrentAxes');
		currentSlice = getappdata(hAxis, 'currentSlice');
		currentTime  = getappdata(hAxis, 'currentTime');

		if (currentTime == 0)
			return
		end
		
		% Create ROI figure if needed; if no ROIs, then quit
		if isempty(roiPlot.hFig) || ~ishandle(roiPlot.hFig)
			bNoRoiData = true;
			for jAxis = reshape(hMainAxes(hMainAxes ~= 0), 1, [])
				roiData = getappdata(jAxis, 'roiData');
				if (numel(roiData{currentSlice}) >= currentTime) && ~isempty(roiData{currentSlice}{currentTime})
					bNoRoiData = false;
				end
			end
			if bNoRoiData
				return
			end
			CreateRoiFigure(hFig);
			roiPlot    = getappdata(hFig, 'roiPlot');
			hTableMean = getappdata(hFig, 'hTableMean');
			hTableStd  = getappdata(hFig, 'hTableStd');
		end

		% Figure out which tab is displayed
		hTabGroup = get(roiPlot.hFig, 'Children');
% 		strSelectedTabTitle = hTabGroup.SelectedTab.Title;
		if ~exist('strSelectedTabTitle', 'var') || isempty(strSelectedTabTitle)
			strSelectedTabTitle = get(get(hTabGroup, 'SelectedTab'), 'Title');
		end

		% Update each subplot
		for iRow = 1:size(roiPlot.hAxes,1)
			for iCol = 1:size(roiPlot.hAxes,2)
				if hMainAxes(iRow,iCol) == 0
					continue
				end

				roi     = getappdata(hMainAxes(iRow,iCol), 'roi');

				switch strSelectedTabTitle
					case 'Mean'
						roiData = getappdata(hMainAxes(iRow,iCol), 'roiData');

						% Remove table data
						if isempty(roi{currentSlice}{currentTime}) && ~isempty(set(hTableMean(iRow,iCol), 'Data'))
							set(hTableMean(iRow,iCol), 'Data', {[]});
							continue
						end

						% Convert roiData to a flat 2D cell array
						aRois = zeros(size(roiData{currentSlice}));
						for iTime = 1:numel(roiData{currentSlice})
							if isempty(roiData{currentSlice}{iTime}) || all(cellfun(@(x) isempty(x), roiData{currentSlice}{iTime}))
								continue
							end
							aRois(iTime) = find(cellfun(@(x) ~isempty(x), roiData{currentSlice}{iTime}), 1, 'last');
						end
						maxRois = max(aRois);
						maxTimes = numel(roiData{currentSlice});

						tableData = cell(maxRois, maxTimes);

						for iTime = 1:maxTimes
							for iRoi = 1:maxRois
								if iRoi <= numel(roiData{currentSlice}{iTime})
									tableData{iRoi,iTime} = roiData{currentSlice}{iTime}{iRoi};
								end
							end
						end

						colNames  = sprintfc('Time %d', 1:maxTimes);
						rowNames  = sprintfc('ROI %d', 1:maxRois);
						colFormat = repmat({'shortg'}, [1 maxTimes]);

						% Update data table
						set(hTableMean(iRow,iCol), 'Data',         tableData)
						set(hTableMean(iRow,iCol), 'ColumnName',   colNames)
						set(hTableMean(iRow,iCol), 'RowName',      rowNames)
						set(hTableMean(iRow,iCol), 'ColumnFormat', colFormat)

					case 'Std Dev'
						roiStd  = getappdata(hMainAxes(iRow,iCol), 'roiStd');

						% Remove table data
						if isempty(roi{currentSlice}{currentTime}) && ~isempty(set(hTableStd(iRow,iCol), 'Data'))
							set(hTableStd(iRow,iCol), 'Data', {[]});
							continue
						end

						% Convert roiStd to a flat 2D cell array
						aRois = zeros(size(roiStd{currentSlice}));
						for iTime = 1:numel(roiStd{currentSlice})
							if isempty(roiStd{currentSlice}{iTime}) || all(cellfun(@(x) isempty(x), roiStd{currentSlice}{iTime}))
								continue
							end
							aRois(iTime) = find(cellfun(@(x) ~isempty(x), roiStd{currentSlice}{iTime}), 1, 'last');
						end
						maxRois = max(aRois);
						maxTimes = numel(roiStd{currentSlice});

						tableData = cell(maxRois, maxTimes);

						for iTime = 1:maxTimes
							for iRoi = 1:maxRois
								if iRoi <= numel(roiStd{currentSlice}{iTime})
									tableData{iRoi,iTime} = roiStd{currentSlice}{iTime}{iRoi};
								end
							end
						end

						colNames  = sprintfc('Time %d', 1:maxTimes);
						rowNames  = sprintfc('ROI %d', 1:maxRois);
						colFormat = repmat({'shortg'}, [1 maxTimes]);

						% Update data table
						set(hTableStd(iRow,iCol), 'Data',         tableData)
						set(hTableStd(iRow,iCol), 'ColumnName',   colNames)
						set(hTableStd(iRow,iCol), 'RowName',      rowNames)
						set(hTableStd(iRow,iCol), 'ColumnFormat', colFormat)
				end
			end
		end
	end
%% -----------------------------------------------------------------------------
function InterpolateRoi(CurrentAxes, aRoi, bUpdateFigure)
% Interpolate an ROI to have a specified number of vertices

	hFig         = get(CurrentAxes, 'Parent');
	prefs        = getappdata(hFig, 'prefs');
	staticRoi    = getappdata(hFig, 'staticRoi');

	currentSlice = getappdata(CurrentAxes, 'currentSlice');
	currentTime  = getappdata(CurrentAxes, 'currentTime');
	roi          = getappdata(CurrentAxes, 'roi');
	
	tol          = 1e-6; % distances smaller than this are considered negligble

	nIgnoredSegments = prefs.interpIgnoreLast;

	if ~exist('aRoi', 'var') || isempty(aRoi)
		aRoi = 1:numel(roi{currentSlice}{currentTime});
	end

	for jRoi = aRoi
		% ROI doesn't exist
		if (jRoi > numel(roi{currentSlice}{currentTime})) || isempty(roi{currentSlice}{currentTime}{jRoi})
			continue
		end

		tmpRoi = roi{currentSlice}{currentTime}{jRoi};

		% A NaN vertex indicates a multi-segment line -- remove it and don't close
		if all(isnan(tmpRoi(end,:)))
			tmpRoi = tmpRoi(1:end-1,:);
			bIsLine = true;
		else
			bIsLine = false;
			% Close the polygon if necessary
			if (sum(abs(diff(tmpRoi([1 end],:)))) > tol)
				tmpRoi = tmpRoi([1:end 1],:);
			end
		end
		
		% Remove duplicate points
		indsBad = sum(abs(diff(tmpRoi)),2) < tol;
		tmpRoi(indsBad,:) = [];

		% Special exception: If we have a line (2 points), and nIgnoredSegments is
		% 0, actually unclose the polygon, otherwise we have two interpolated lines
		% on top of each other
		if (size(tmpRoi,1) == 3) && (nIgnoredSegments == 0)
			tmpRoi = tmpRoi(1:end-1,:);
		else
			tmpRoi = tmpRoi(1:end-nIgnoredSegments,:);
		end

		% Insufficient number of points to interpolate (e.g. a point)
		if size(tmpRoi,1) < 2
			continue
		end
		
		% Interpolate x and y independently, based on cumulative length along the polygon's circumference
		diffDist = sqrt(diff(tmpRoi(:,1)).^2 + diff(tmpRoi(:,2)).^2);
		cumDist = cat(1, 0, cumsum(diffDist));

		% Number of vertices to interpolate to
		nPts = prefs.interpNPts + 1 - nIgnoredSegments; % We add one at the end to complete the polygon

		% Number of points necessary to have distance between points equal to the
		% minimum distance between any two existing vertices.
		nMinPts = ceil(cumDist(end) / min(diffDist));

		nPts = max(nPts, nMinPts);

		% This could potentially lead to a very large number of points if there's
		% even one complex shape, so cap it
		nPts = min(nPts, 50);		
		
		tmpInterpRoi = zeros(nPts,2);
		tmpInterpRoi(:,1) = interp1(cumDist, tmpRoi(:,1), linspace(0,cumDist(end),nPts), prefs.interpAlgorithm);
		tmpInterpRoi(:,2) = interp1(cumDist, tmpRoi(:,2), linspace(0,cumDist(end),nPts), prefs.interpAlgorithm);

		if ~bIsLine && (size(tmpRoi,1) ~= 2) && (nIgnoredSegments == 0)
			% Remove vertex used to close polygon
			tmpInterpRoi = tmpInterpRoi(1:end-1,:);
		else
			tmpInterpRoi = cat(1, tmpInterpRoi, roi{currentSlice}{currentTime}{jRoi}(end-nIgnoredSegments+2:end,:));
		end

		% Re-add NaN marker for line segments
		if bIsLine
			tmpInterpRoi = cat(1, tmpInterpRoi, [nan nan]);
		end

		roi{currentSlice}{currentTime}{jRoi} = tmpInterpRoi;
	end

	setappdata(CurrentAxes, 'roi', roi);
	DoStaticAndTranscribeAnd3D(hFig, aRoi);
	
% 		if staticRoi
% 			for iTime = 1:numel(roi{currentSlice})
% 				roi{currentSlice}{iTime}{jRoi} = tmpInterpRoi;
% 			end
% 		else
% 			roi{currentSlice}{currentTime}{jRoi} = tmpInterpRoi;
% 		end
% 	end
% 
% 	if ~exist('bUpdateFigure', 'var') || bUpdateFigure
% 		for jRoi = aRoi
% 			UpdateRoiData(hFig, jRoi);
% 		end
% 	end
end

% % function InterpolateRoi(hFig)
% % % Interpolate an ROI to have a specified number of vertices
% % 	nIgnoredSegments = 2;
% % 	nPts = 50;
% % 	nPts = nPts + 1; % We add one at the end to complete the polygon
% % 
% % 	CurrentAxes  = get(hFig, 'CurrentAxes');
% % 	selectedRoi  = getappdata(hFig, 'selectedRoi');
% % 	currentSlice = getappdata(CurrentAxes, 'currentSlice');
% % 	currentTime  = getappdata(CurrentAxes, 'currentTime');
% % 	staticRoi    = getappdata(hFig, 'staticRoi');
% % 
% % % 	if isempty(selectedRoi)
% % % 		return
% % % 	end
% % 
% % 	roi = getappdata(CurrentAxes, 'roi');
% % 
% % 	% This is necessary because we allow an ROI number to be selected, even if
% % 	% it hasn't yet been drawn
% % % 	if selectedRoi > numel(roi{currentSlice}{currentTime})
% % % 		return
% % % 	end
% % 
% % 	for iRoi = 1%selectedRoi%1:numel(roi{currentSlice}{currentTime})
% % 		tmpRoi = roi{currentSlice}{currentTime}{iRoi};
% % 		tmpRoi = cat(1, tmpRoi, tmpRoi(1,:));		
% % 		tmpRoi = tmpRoi(1:end-nIgnoredSegments,:);
% % 
% % 		diffDist = sqrt(diff(tmpRoi(:,1)).^2 + diff(tmpRoi(:,2)).^2);
% % 		cumDist = cat(1, 0, cumsum(diffDist));
% % 
% % 		tmpInterpRoi = zeros(nPts,2);
% % 		tmpInterpRoi(:,1) = interp1(cumDist, tmpRoi(:,1), linspace(0,cumDist(end),nPts));
% % 		tmpInterpRoi(:,2) = interp1(cumDist, tmpRoi(:,2), linspace(0,cumDist(end),nPts));
% % 
% % 		if (nIgnoredSegments == 0)
% % 			tmpInterpRoi = tmpInterpRoi(1:end-1,:);
% % 		else
% % 			tmpInterpRoi = cat(1, tmpInterpRoi, roi{currentSlice}{currentTime}{iRoi}(end-nIgnoredSegments+2:end,:));
% % 		end
% % 		if staticRoi
% % 			for iTime = 1:numel(roi{currentSlice})
% % 				roi{currentSlice}{iTime}{iRoi} = tmpInterpRoi;
% % 			end
% % 		else
% % 			roi{currentSlice}{currentTime}{iRoi} = tmpInterpRoi;
% % 		end
% % 	end
% % 
% % 	setappdata(CurrentAxes, 'roi', roi);
% % % 	UpdateRoiData(hFig, selectedRoi);
% % 	UpdateRoiData(hFig, 1);
% % end

%% -----------------------------------------------------------------------------
function InterpolateAcrossTime(hFig)
	CurrentAxes = get(hFig, 'CurrentAxes');
	fcnHandles = getappdata(hFig, 'fcnHandles');

	roi = getappdata(CurrentAxes, 'roi');

	aNonEmptyTimes = find(cellfun(@(x) ~isempty(x), roi{1}));
	aNonEmptyTimes(cellfun(@(x) all(cellfun(@(y) isempty(y), x)), roi{1}(aNonEmptyTimes))) = [];  % Ignore times with blank ROIs

	interpRoi = {cell(size(roi{1}))};
	for iRoi = 1:max(cellfun(@(x) numel(x), roi{1}(aNonEmptyTimes)))
		% Ignore times that don't have the right number of ROIs
		aGoodTimes = aNonEmptyTimes(cellfun(@(x) numel(x), roi{1}(aNonEmptyTimes)) >= iRoi);

		% Ignore times with empty ROIs
		aGoodTimes = aGoodTimes(cellfun(@(x) ~isempty(x{iRoi}), roi{1}(aGoodTimes)));

		% Ignore times with a placeholder ROI, but make sure to put that ROI back into interpRoi
		aPlaceholderInds = find(cellfun(@(x) numel(x{iRoi}) == 2, roi{1}(aGoodTimes)));
		for jInd = aPlaceholderInds
			interpRoi{1}{aGoodTimes(jInd)}{iRoi} = roi{1}{aGoodTimes(jInd)}{iRoi};
		end
		aGoodTimes(aPlaceholderInds) = [];

		% Nothing to do!
		if isempty(aGoodTimes)
			continue
		end

		% Ensure that every ROI has the same number of vertices
		cGoodRois = permute(cellfun(@(x) x{iRoi}, roi{1}(aGoodTimes), 'UniformOutput', false), [1 3 2]);
		if (numel(unique(cellfun(@(x) size(x,1), cGoodRois))) ~= 1)
			% Can't interpolate, but make sure the existing ROIs are placed back into interpRoi
			for jInd = aGoodTimes
				interpRoi{1}{jInd}{iRoi} = roi{1}{jInd}{iRoi};
			end
			continue
		end

		% Extract the existing ROIs into a 3D matrix
		goodRois = cell2mat(cGoodRois);

		% These are the times we need interpolated ROIs for
		aInterpTimes = aGoodTimes(1):aGoodTimes(end);
		allRois = zeros([size(goodRois,1) 2 numel(aInterpTimes)]);

		% Interpolate!
		for iVert = 1:size(goodRois,1)
			allRois(iVert,1,:) = interp1(aGoodTimes, squeeze(goodRois(iVert,1,:))', aInterpTimes);
			allRois(iVert,2,:) = interp1(aGoodTimes, squeeze(goodRois(iVert,2,:))', aInterpTimes);
		end

		% Re-format back for SimpleViewer
		for iTime = 1:numel(aInterpTimes)
			interpRoi{1}{aInterpTimes(iTime)}{iRoi} = allRois(:,:,iTime);
		end
	end

	% Save the new ROIs back to figure and update the display
	setappdata(CurrentAxes, 'roi', interpRoi);
	fcnHandles.UpdateAllRoiHandles(hFig);
end

%% -----------------------------------------------------------------------------
function StartPanZoom(hFig, evt)
% Called when a mouse-click is started; dispatches to the appropriate pan/zoom function

	% Update the selected subplot, ROI
	MouseDownSelect(hFig, evt);

	% If the click occurs on an existing ROI object, abort and call that object's callback
	hObj  = get(hFig, 'CurrentObject');
	hFcn  = get(hObj, 'ButtonDownFcn');
	if ~isempty(hFcn)
		MouseDownSelect(hFig, evt);
		hFcn(hFig, evt);
		return
	end

	hAxis = get(hFig, 'CurrentAxes');
	switch(get(hFig, 'SelectionType'))
		case 'normal'
			% Left Click: Pan
			setappdata(hFig, 'InitPoint', get(hAxis, 'CurrentPoint'))
			set(hFig, 'WindowButtonMotionFcn', @AdjustPan);
			set(hFig, 'WindowButtonUpFcn',     @EndPanZoom);
			set(hFig, 'Pointer', 'fleur')
		case 'alt'
			% Right Click: Zoom
% 		szFig = get(hFig, 'Position');
% 		posSubplot = get(hAxis, 'Position') .* szFig;
% 		XLim = mean(get(iAxes, 'XLim')) + [-posSubplot(3)/newZoom posSubplot(3)/newZoom] ./ 2;
% 		YLim = mean(get(iAxes, 'YLim')) + [-posSubplot(4)/newZoom posSubplot(4)/newZoom] ./ 2;

			szImg = size(get(findobj(get(hAxis, 'Children'), 'Type', 'image'), 'CData'));
			szImg = szImg(1:2);
			szAxis = [diff(get(hAxis, 'YLim')) diff(get(hAxis, 'XLim'))];
			currentZoom = max(szImg ./ szAxis);

			setappdata(hFig, 'InitPoint', get(hAxis, 'CurrentPoint'))
			setappdata(hFig, 'InitZoom',  currentZoom)
			set(hFig, 'WindowButtonMotionFcn', @AdjustZoom);
			set(hFig, 'WindowButtonUpFcn',     @EndPanZoom);
			setptr(hFig, 'glass')
		case 'open'
			% Double Click: Reset to auto window level for each subplot
			ResetPanZoom(hFig);
	end
end

% ------------------------------------------------------------------------------
function EndPanZoom(hFig, evt)
% Called when a mouse-click is lifted; stops mouse movements from manually adjusting pan

	pData = getappdata(hFig, 'pData');
	set(hFig, 'Pointer', 'custom', 'PointerShapeCData', pData.panZoom, 'PointerShapeHotSpot', [8 8])

	set(hFig, 'WindowButtonMotionFcn', '');
	set(hFig, 'WindowButtonUpFcn',     '');
	
	% New auto window level mode
	bAutoWindowLevel = getappdata(hFig, 'bAutoWindowLevel');
	if ~isempty(bAutoWindowLevel) && bAutoWindowLevel
		AutoWindowLevel(hFig);
	end
end

% ------------------------------------------------------------
function AdjustZoom(hFig, evt, zoomFactor)
% Mouse up/down adjusts zoom (for all axes)

	CurrentAxes = get(hFig, 'CurrentAxes');
	CurrentModifier = get(hFig, 'CurrentModifier');

	InitPoint = getappdata(hFig, 'InitPoint');
	InitZoom  = getappdata(hFig, 'InitZoom');

	% Abort the figure has never been zoomed before (happens when the ResizeFcn is called)
	if isempty(InitPoint) || isempty(InitZoom)
		return
	end
	
	CurrentPoint = get(CurrentAxes, 'CurrentPoint');
	szImg = size(get(findobj(get(CurrentAxes, 'Children'), 'Type', 'image'), 'CData'));

	% Use fraction  i.e. relative to position to the originally clicked point
	% to determine the change in window and level
	delta = CurrentPoint(1,1:2) - InitPoint(1,1:2);

% 	[InitZoom delta(2)/szImg(1)]
	
	newZoom = InitZoom / 1.1^(15 * delta(2)/szImg(1));
	
	if exist('zoomFactor', 'var')
		newZoom = InitZoom * zoomFactor;
	end

	% Set some limits to zoom
	if newZoom < 1
		newZoom = 1;
	elseif newZoom > 10
		newZoom = 10;
	end
	
%	newZoom = 2;
	
	szFig = get(hFig, 'Position');

	% Set for all axes
	hAxes = getappdata(hFig, 'hAxes');
	hAxes = reshape(hAxes(hAxes ~= 0), 1, []);
	for iAxes = hAxes
		% 'Alt' allows pan/zoom to affect only one subplot (probably a cleaner way to implement this though)
		if (~isempty(CurrentModifier) && strcmp(CurrentModifier, 'alt')) && (iAxes ~= CurrentAxes)
			continue
		end

		szImg = size(get(findobj(get(iAxes, 'Children'), 'Type', 'image'), 'CData'));
		posSubplot = get(iAxes, 'Position') .* szFig;

		% Don't adjust for empty subplots
		if (all(szImg) == 0)
			continue
		end
		
		currXLim = get(iAxes, 'XLim') - 0.5;
		currYLim = get(iAxes, 'YLim') - 0.5;

		if szImg(1) >= szImg(2)
% % 			XLim = szImg(1)/2 + [-1 1]*(szImg(1)/2)/newZoom + 0.5;
% % 			YLim = szImg(2)/2 + [-1 1]*(szImg(1)*posSubplot(2)/posSubplot(3)/2)/newZoom + 0.5;
% 			YLim = szImg(1)/2 + [-1 1]*(szImg(1)/2)/newZoom + 0.5;
% 			XLim = szImg(2)/2 + [-1 1]*(szImg(1)*posSubplot(3)/posSubplot(4)/2)/newZoom + 0.5;
			YLim = mean(currYLim) + [-1 1]*(szImg(1)/2)/newZoom + 0.5;
			XLim = mean(currXLim) + [-1 1]*(szImg(1)*posSubplot(3)/posSubplot(4)/2)/newZoom + 0.5;
			
			XLim(1) = max([XLim(1) 0.5]);
			XLim(2) = min([XLim(2) szImg(2)+0.5]);
		elseif szImg(2) > szImg(1)
% % 			YLim = szImg(1)/2 + [-1 1]*(szImg(1)/2)/newZoom + 0.5;
% % 			XLim = szImg(2)/2 + [-1 1]*(szImg(1)*posSubplot(4)/posSubplot(3)/2)/newZoom + 0.5;
% 			XLim = szImg(2)/2 + [-1 1]*(szImg(2)/2)/newZoom + 0.5;
% 			YLim = szImg(1)/2 + [-1 1]*(szImg(2)*posSubplot(4)/posSubplot(3)/2)/newZoom + 0.5;
			XLim = mean(currXLim) + [-1 1]*(szImg(2)/2)/newZoom + 0.5;
			YLim = mean(currYLim) + [-1 1]*(szImg(2)*posSubplot(4)/posSubplot(3)/2)/newZoom + 0.5;

			YLim(1) = max([YLim(1) 0.5]);
			YLim(2) = min([YLim(2) szImg(1)+0.5]);
		end

		
% 		XLim = mean(get(iAxes, 'XLim')) + [-szImg(2)/newZoom szImg(2)/newZoom] ./ 2;
% 		YLim = mean(get(iAxes, 'YLim')) + [-szImg(1)/newZoom szImg(1)/newZoom] ./ 2;

% 		XLim = mean(get(iAxes, 'XLim')) + [-posSubplot(3)/newZoom posSubplot(3)/newZoom] ./ 2;
% 		YLim = mean(get(iAxes, 'YLim')) + [-posSubplot(4)/newZoom posSubplot(4)/newZoom] ./ 2;

		set(iAxes, 'XLim', XLim);
		set(iAxes, 'YLim', YLim);
	end

	% Save the new zoom state, as it's used to re-zoom when the window size is adjusted
% 	setappdata(hFig, 'InitZoom', newZoom);
	
	% When true, SetTimeSlice won't update XLim, YLim
	setappdata(hFig, 'bManualZoom', true);
end

% ------------------------------------------------------------
function AdjustPan(hFig, evt)
% Left/right adjusts window, up/down adjusts level
% Adapted from 'imagescn'.

	CurrentAxes = get(hFig, 'CurrentAxes');
	CurrentModifier = get(hFig, 'CurrentModifier');

	InitPoint = getappdata(hFig, 'InitPoint');
	CurrentPoint = get(CurrentAxes, 'CurrentPoint');

	XLim = get(CurrentAxes, 'XLim');
	YLim = get(CurrentAxes, 'YLim');
	szImg = size(get(findobj(get(CurrentAxes, 'Children'), 'Type', 'image'), 'CData'));

	% Use fraction  i.e. relative to position to the originally clicked point
	% to determine the change in window and level
	delta = CurrentPoint(1,1:2) - InitPoint(1,1:2);

	if (XLim(1) - delta(1)) < 0.5
		delta(1) = 0.5 + XLim(1);
	elseif (XLim(2) - delta(1)) > szImg(2) + 0.5
		delta(1) = XLim(2) - (szImg(2) + 0.5);
	end

	if (YLim(1) - delta(2)) < 0.5
		delta(2) = 0.5 + YLim(1);
	elseif (YLim(2) - delta(2)) > szImg(1) + 0.5
		delta(2) = YLim(2) - (szImg(1) + 0.5);
	end

	% Set for all axes
	hAxes = getappdata(hFig, 'hAxes');
	for iAxes = reshape(hAxes(hAxes ~= 0), 1, [])
		% 'Alt' allows pan/zoom to affect only one subplot (probably a cleaner way to implement this though)
		if (~isempty(CurrentModifier) && strcmp(CurrentModifier, 'alt')) && (iAxes ~= CurrentAxes)
			continue
		end

		XLim = get(iAxes, 'XLim');
		YLim = get(iAxes, 'YLim');
			set(iAxes, 'XLim', XLim - delta(1));
		set(iAxes, 'YLim', YLim - delta(2));
	end
	% When true, SetTimeSlice won't update XLim, YLim
	setappdata(hFig, 'bManualZoom', true);
end


% ------------------------------------------------------------------------------
function ResetPanZoom(hFig, evt)
% Each subplot is independently automatically window leveled
	hAxes = getappdata(hFig, 'hAxes');
	
	for hAxis = reshape(hAxes(hAxes ~= 0), 1, [])
		szImg  = size(get(getappdata(hAxis, 'hImage'), 'CData'));
		if all(szImg == 0)
			continue
		end
		set(hAxis, 'XLim', [0.5 szImg(2)+0.5])
		set(hAxis, 'YLim', [0.5 szImg(1)+0.5])
	end

	% When false, SetTimeSlice update XLim, YLim for each subplot automatically
	% (and independently)
	setappdata(hFig, 'bManualZoom', false);
end













%% -----------------------------------------------------------------------------
function ToggleRoiDisplay(hFig, evt)
% Turns all ROI handles on/off

	hideRois = getappdata(hFig, 'hideRois');

	if hideRois
		strVisible = 'on';
	else
		strVisible = 'off';
	end

	hAxes = getappdata(hFig, 'hAxes');

	for hAxis = hAxes(:)'
		roiHandles = getappdata(hAxis, 'roiHandles');
		for iRoi = 1:numel(roiHandles)
			hAll = [roiHandles(iRoi).line roiHandles(iRoi).vertices roiHandles(iRoi).center roiHandles(iRoi).label];
			set(hAll, 'Visible', strVisible)
		end
	end

	setappdata(hFig, 'hideRois', ~hideRois);

	% Apply to SimpleViewer3D, if applicable
	hSV3D = getappdata(hFig, 'hSV3D');
	if ~isempty(hSV3D) && ishandle(hSV3D)
		fcnHandles = getappdata(hSV3D, 'fcnHandles');
		fcnHandles.ToggleRoiDisplay(hSV3D)
	end
end

%% -----------------------------------------------------------------------------
function GotoHomeTime(hFig)
% Go to home time frames in order
	currentAxes  = get(hFig,               'CurrentAxes');
	prefs        = getappdata(hFig,        'prefs');
	currentSlice = getappdata(currentAxes, 'currentSlice');
	currentTime  = getappdata(currentAxes, 'currentTime');

	iHomeTime = find(prefs.homeTime == currentTime);

	if isempty(iHomeTime)
		iHomeTime = 0;
	end

	newTime = mod(iHomeTime+1-1, numel(prefs.homeTime))+1;
	SetTimeSlice(hFig, prefs.homeTime(newTime), currentSlice);
end

%% -----------------------------------------------------------------------------
function SetTimeSlice(hFig, newTime, newSlice, bForceReload, bUpdateSingleSubplot)
% Sets the time/slice for all subplots within a figure
% 	if isMultipleCall(); return; end
	prefs       = getappdata(hFig,        'prefs');
	persistent lastCallTime
	if ~isempty(lastCallTime)
		if etime(clock, lastCallTime) < 1/prefs.maxFrameRate
% 			disp('Exeeding frame rate limiter!')
			return
		end
	end

	currentAxes = get(hFig,               'CurrentAxes');
	hAxes       = getappdata(hFig,        'hAxes');
% 	prefs       = getappdata(hFig,        'prefs');
	bManualZoom = getappdata(hFig,        'bManualZoom');
	imgData     = getappdata(currentAxes, 'imgData');
	oldTime     = getappdata(currentAxes, 'currentTime');
	oldSlice    = getappdata(currentAxes, 'currentSlice');

	if isempty(imgData)
		return
	end

	if (newTime == oldTime) && (newSlice == oldSlice) && ~(exist('bForceReload', 'var') && ~isempty(bForceReload) && bForceReload)
		return
	end

	% Allow only a single subplot to be updated
	if ~exist('bUpdateSingleSubplot', 'var') || isempty(bUpdateSingleSubplot)
		bUpdateSingleSubplot = false;
	end
	
	% ROI drawing places some limitations on what you can change
	if ~isempty(getappdata(currentAxes, 'inProgressRoi'))
		if newSlice ~= oldSlice
			set(getappdata(hFig, 'hStatusText'), 'String', 'Error: Cannot change slice while drawing an ROI', 'Color', 'r');
			return
		end
		if ~getappdata(hFig, 'staticRoi') && (newTime ~= oldTime)
			set(getappdata(hFig, 'hStatusText'), 'String', 'Error: Cannot change time while drawing an ROI if staticRoi is off', 'Color', 'r');
			return
		end
	end
	
	% Auto-resize means that the window size will change to keep the same zoom ratio as the previous frame
	if prefs.autoResize
		scale = CalculateCurrentScale(hFig);
		if scale < 0.5
			scale = 0.5;  % Don't allow the window to become too small
		end
	end

	if (newSlice ~= oldSlice)
		% When changing slice, "clip" the time if necessary
		newSlice = mod(newSlice-1, length(imgData))+1;
		newTime  = min([length(imgData{newSlice}) newTime]);
	elseif (newTime ~= oldTime)
 		% Slice is already valid
		newTime  = mod(newTime -1, length(imgData{newSlice}))+1;
	end

	% Update all subplots
	for jAxis = reshape(hAxes(hAxes ~= 0), 1, [])
		currSubplot   = find(hAxes == jAxis, 1);

		if bUpdateSingleSubplot
			if (jAxis ~= currentAxes)
				continue
			end
		end

		newSubTime    = newTime;
		newSubSlice   = newSlice;
		imgData       = getappdata(jAxis, 'imgData');
		txtData       = getappdata(jAxis, 'txtData');
		roi           = getappdata(jAxis, 'roi');
		roiHandles    = getappdata(jAxis, 'roiHandles');
		hImage        = getappdata(jAxis, 'hImage');
		hText         = getappdata(jAxis, 'hText');
		inProgressRoi = getappdata(jAxis, 'inProgressRoi');

		if isempty(imgData)
			continue
		end

		% Range check (these are clipped, not cyclic)
		newSubSlice = min([newSubSlice length(imgData             )]);
		newSubTime  = min([newSubTime  length(imgData{newSubSlice})]);

		if (newSubSlice == getappdata(jAxis, 'currentSlice')) && (newSubTime  == getappdata(jAxis, 'currentTime'))
			continue
		end

		% Image
		oldImgSize = size(get(hImage, 'CData'));
		newImgSize = size(imgData{newSubSlice}{newSubTime});
		set(hImage, 'CData',  imgData{newSubSlice}{newSubTime});

		if isempty(bManualZoom) || ~bManualZoom
			if ~isempty(imgData{newSubSlice}{newSubTime}) && ((oldImgSize(1) ~= newImgSize(1)) || (oldImgSize(2) ~= newImgSize(2)))
				set(jAxis, 'XLim', [0.5 newImgSize(2)+0.5])
				set(jAxis, 'YLim', [0.5 newImgSize(1)+0.5])
			end
		end

		% Text
		strText = ComposeText(txtData{newSubSlice}{min([newSubTime  length(txtData{newSubSlice})])}, ...
													newSubSlice, length(imgData), ...
													newSubTime,  length(imgData{newSubSlice}));
		set(hText, 'String', strText);

		% This has to be done before calling UpdateRoiHandles()
		setappdata(jAxis, 'currentSlice', newSubSlice);
		setappdata(jAxis, 'currentTime',  newSubTime );
		setappdata(jAxis, 'currentSliceImg', []);
		setappdata(jAxis, 'currentTimeImg',  []);

		if ~isempty(inProgressRoi)
			continue
		end

		% Update ROI handles for rois on this time/slice
		for iRoi = 1:length(roi{newSubSlice}{newSubTime})
% 			% TODO: this makes it slower,  but need this to update roiPerimeter/roiArea
% 			UpdateRoiData(hFig, iRoi);
			UpdateRoiHandles(jAxis, iRoi)
		end

		% Hide any additional plot handles
		for iRoi = length(roi{newSubSlice}{newSubTime})+1:length(roiHandles)
			set(roiHandles(iRoi).line,     'XData', [],   'YData', []);
			set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
			set(roiHandles(iRoi).center,   'XData', [],   'YData', []);
			if ~strcmp(prefs.roiLabelMode, 'none')
				set(roiHandles(iRoi).label,    'Position', nan(1,2));
			end
		end

		if prefs.useDicomColormaps
			% Use custom colour map if present
			extras = getappdata(hFig, 'extras');
			if isfield(extras, 'headers')
				info = extras.headers{currSubplot}{newSubSlice}{newSubTime};
				if all(isfield(info, {'RedPaletteColorLookupTableData', 'GreenPaletteColorLookupTableData', 'BluePaletteColorLookupTableData', 'RedPaletteColorLookupTableDescriptor'}))
					cmap = cat(2, info.RedPaletteColorLookupTableData, info.GreenPaletteColorLookupTableData, info.BluePaletteColorLookupTableData);
					if isrow(cmap)
						cmap = reshape(cmap, [], 3);
					end
					cmap = double(cmap) / double(max(cmap(:)));

					% Assume this is the same as Green and Blue's palette
					% The first value is the number of entries in the lookup table. The second value is the first input value mapped.
					CLim = double(info.RedPaletteColorLookupTableDescriptor(2) + [0 info.RedPaletteColorLookupTableDescriptor(1)-1]);

					% But.. MATLAB's 'painters' renderer doesn't support more than 256 colormap entries
					t = linspace(CLim(1), CLim(2), 256)';
					cmap = cat(2, interp1(CLim(1):CLim(2), cmap(:,1), t), ...
												interp1(CLim(1):CLim(2), cmap(:,2), t), ...
												interp1(CLim(1):CLim(2), cmap(:,3), t));
				else
					cmap = gray(256);
				end
				% Preset window level
				if all(isfield(info, {'WindowCenter', 'WindowWidth'}))
					CLim = info.WindowCenter(1) + info.WindowWidth/2*[-1 1];
					set(jAxis, 'CLim', CLim)
				end
				colormap(jAxis, cmap)
			end
		end
	end


	setappdata(currentAxes, 'lastTime',  oldTime);
	setappdata(currentAxes, 'lastSlice', oldSlice);

	% New auto window level mode
	bAutoWindowLevel = getappdata(hFig, 'bAutoWindowLevel');
	if ~isempty(bAutoWindowLevel) && bAutoWindowLevel
		AutoWindowLevel(hFig);
	end
	
	% FIXME: This should be a remnant of when drawing points and it required changing the current axis
% 	set(hFig, 'CurrentAxes', currentAxes);
	if prefs.autoResize
		SetFigureScale(hFig, scale);
	end
% 	drawnow expose
% 	pause(1)

% 	if (newSlice ~= oldSlice)
		UpdateRoiFigure(hFig);
% 	end

	% Update intersection plots
	UpdatePlotIntersect(hFig);

	% Apply to SimpleViewer3D, if applicable
	hSV3D    = getappdata(hFig, 'hSV3D');
	hSV3DSet = getappdata(hFig, 'hSV3DSet');	

	if isempty(hSV3DSet)
		hSV3DSet = 1;
	end

	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		% First, set it to the first "set" of images
% 		hSurf = getappdata(hSV3D, 'hSurf');
% 		setappdata(hSV3D, 'lastCurrentObj', hSurf(1));

		fcnHandles = getappdata(hSV3D, 'fcnHandles');
		fcnHandles.SetTimeSlice(hSV3D, newTime, newSlice, hSV3DSet)
	end

	% Apply to another SimpleViewer, if applicable
	hSV = getappdata(hFig, 'hSV');
	if ~isempty(hSV) && ishandle(hSV)
		fcnHandles = getappdata(hSV, 'fcnHandles');
		fcnHandles.SetTimeSlice(hSV, newTime, newSlice)
	end

	lastCallTime = clock;
% 	pause(1/30)
end

%% -----------------------------------------------------------------------------
function SetTimeSliceStruct(hFig, newTime, newSlice)
% Sets the time/slice for all subplots within a figure

% For some reason this is much slower than SetTimeSlice if only the current
% frame is being updated and the others are not.  FIXME

% 	if isMultipleCall(); return; end
	prefs       = getappdata(hFig,        'prefs');
	persistent lastCallTime
	if ~isempty(lastCallTime)
		if etime(clock, lastCallTime) < 1/prefs.maxFrameRate
% 			disp('Exeeding frame rate limiter!')
			return
		end
	end
	
	if verLessThan('matlab','8.4.0')
		SetTimeSlice(hFig, newTime, newSlice);
		return
	end
	
	currentAxes = get(hFig,               'CurrentAxes');
	hAxes       = getappdata(hFig,        'hAxes');
% 	prefs       = getappdata(hFig,        'prefs');
	bManualZoom = getappdata(hFig,        'bManualZoom');
	imgData     = getappdata(currentAxes, 'imgData');
	oldTime     = getappdata(currentAxes, 'currentTime');
	oldSlice    = getappdata(currentAxes, 'currentSlice');

	if isempty(imgData)
		return
	end

	if (newTime == oldTime) && (newSlice == oldSlice)
		return
	end

	% ROI drawing places some limitations on what you can change
	if ~isempty(getappdata(currentAxes, 'inProgressRoi'))
		if newSlice ~= oldSlice
			set(getappdata(hFig, 'hStatusText'), 'String', 'Error: Cannot change slice while drawing an ROI', 'Color', 'r');
			return
		end
		if ~getappdata(hFig, 'staticRoi') && (newTime ~= oldTime)
			set(getappdata(hFig, 'hStatusText'), 'String', 'Error: Cannot change time while drawing an ROI if staticRoi is off', 'Color', 'r');
			return
		end
	end
	
	% Auto-resize means that the window size will change to keep the same zoom ratio as the previous frame
	if prefs.autoResize
		scale = CalculateCurrentScale(hFig);
		if scale < 0.5
			scale = 0.5;  % Don't allow the window to become too small
		end
	end

	if (newSlice ~= oldSlice)
		% When changing slice, "clip" the time if necessary
		newSlice = mod(newSlice-1, length(imgData))+1;
		newTime  = min([length(imgData{newSlice}) newTime]);
	elseif (newTime ~= oldTime)
 		% Slice is already valid
		newTime  = mod(newTime -1, length(imgData{newSlice}))+1;
	end

	% Update all subplots
	for jAxis = reshape(hAxes(hAxes ~= 0), 1, [])
		currSubplot   = find(hAxes == jAxis, 1);
		newSubTime    = newTime;
		newSubSlice   = newSlice;
% 		imgData       = getappdata(jAxis, 'imgData');
% 		txtData       = getappdata(jAxis, 'txtData');
% 		roi           = getappdata(jAxis, 'roi');
% 		roiHandles    = getappdata(jAxis, 'roiHandles');
% 		hImage        = getappdata(jAxis, 'hImage');
% 		hText         = getappdata(jAxis, 'hText');
% 		inProgressRoi = getappdata(jAxis, 'inProgressRoi');

		adAxis = getappdata(jAxis);

		if isempty(imgData)
			continue
		end

		% Range check (these are clipped, not cyclic)
		newSubSlice = min([newSubSlice length(adAxis.imgData             )]);
		newSubTime  = min([newSubTime  length(adAxis.imgData{newSubSlice})]);

		if (newSubSlice == adAxis.currentSlice) && (newSubTime == adAxis.currentTime)
			continue
		end

		% Image
		oldImgSize = size(adAxis.hImage.CData);
		newImgSize = size(adAxis.imgData{newSubSlice}{newSubTime});
		adAxis.hImage.CData = adAxis.imgData{newSubSlice}{newSubTime};

		if isempty(bManualZoom) || ~bManualZoom
			if ~isempty(adAxis.imgData{newSubSlice}{newSubTime}) && ((oldImgSize(1) ~= newImgSize(1)) || (oldImgSize(2) ~= newImgSize(2)))
				set(jAxis, 'XLim', [0.5 newImgSize(2)+0.5])
				set(jAxis, 'YLim', [0.5 newImgSize(1)+0.5])
			end
		end

		% Text
		strText = ComposeText(adAxis.txtData{newSubSlice}{min([newSubTime  length(adAxis.txtData{newSubSlice})])}, ...
													newSubSlice, length(adAxis.imgData), ...
													newSubTime,  length(adAxis.imgData{newSubSlice}));
		adAxis.hText.String = strText;

		% This has to be done before calling UpdateRoiHandles()
		setappdata(jAxis, 'currentSlice', newSubSlice);
		setappdata(jAxis, 'currentTime',  newSubTime );
		setappdata(jAxis, 'currentSliceImg', []);
		setappdata(jAxis, 'currentTimeImg',  []);

		if ~isempty(adAxis.inProgressRoi)
			continue
		end

		% Update ROI handles for rois on this time/slice
		for iRoi = 1:length(adAxis.roi{newSubSlice}{newSubTime})
% 			% TODO: this makes it slower,  but need this to update roiPerimeter/roiArea
% 			UpdateRoiData(hFig, iRoi);
			UpdateRoiHandles(jAxis, iRoi)
		end

		% Hide any additional plot handles
		for iRoi = length(adAxis.roi{newSubSlice}{newSubTime})+1:length(adAxis.roiHandles)
			set(adAxis.roiHandles(iRoi).line,     'XData', [],   'YData', []);
			set(adAxis.roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
			set(adAxis.roiHandles(iRoi).center,   'XData', [],   'YData', []);
			if ~strcmp(prefs.roiLabelMode, 'none')
				set(adAxis.roiHandles(iRoi).label,    'Position', nan(1,2));
			end
		end

		if prefs.useDicomColormaps
			% Use custom colour map if present
			extras = getappdata(hFig, 'extras');
			if isfield(extras, 'headers')
				info = extras.headers{currSubplot}{newSubSlice}{newSubTime};
				if all(isfield(info, {'RedPaletteColorLookupTableData', 'GreenPaletteColorLookupTableData', 'BluePaletteColorLookupTableData', 'RedPaletteColorLookupTableDescriptor'}))
					cmap = cat(2, info.RedPaletteColorLookupTableData, info.GreenPaletteColorLookupTableData, info.BluePaletteColorLookupTableData);
					if isrow(cmap)
						cmap = reshape(cmap, [], 3);
					end
					cmap = double(cmap) / double(max(cmap(:)));

					% Assume this is the same as Green and Blue's palette
					% The first value is the number of entries in the lookup table. The second value is the first input value mapped.
					CLim = double(info.RedPaletteColorLookupTableDescriptor(2) + [0 info.RedPaletteColorLookupTableDescriptor(1)-1]);

					% But.. MATLAB's 'painters' renderer doesn't support more than 256 colormap entries
					t = linspace(CLim(1), CLim(2), 256)';
					cmap = cat(2, interp1(CLim(1):CLim(2), cmap(:,1), t), ...
												interp1(CLim(1):CLim(2), cmap(:,2), t), ...
												interp1(CLim(1):CLim(2), cmap(:,3), t));
				else
					cmap = gray(256);
				end
				% Preset window level
				if all(isfield(info, {'WindowCenter', 'WindowWidth'}))
					CLim = info.WindowCenter(1) + info.WindowWidth/2*[-1 1];
					set(jAxis, 'CLim', CLim)
				end
				colormap(jAxis, cmap)
			end
		end
	end

	setappdata(currentAxes, 'lastTime',  oldTime);
	setappdata(currentAxes, 'lastSlice', oldSlice);

	% New auto window level mode
	bAutoWindowLevel = getappdata(hFig, 'bAutoWindowLevel');
	if ~isempty(bAutoWindowLevel) && bAutoWindowLevel
		AutoWindowLevel(hFig);
	end
	
	% FIXME: This should be a remnant of when drawing points and it required changing the current axis
% 	set(hFig, 'CurrentAxes', currentAxes);
	if prefs.autoResize
		SetFigureScale(hFig, scale);
	end
% 	drawnow expose
% 	pause(1)

% 	if (newSlice ~= oldSlice)
	if ~isempty(adAxis.roi{newSubSlice}{newSubTime})
		UpdateRoiFigure(hFig);
	end
% 	end

	% Update intersection plots
	UpdatePlotIntersect(hFig);

	% Apply to SimpleViewer3D, if applicable
	hSV3D    = getappdata(hFig, 'hSV3D');
	hSV3DSet = getappdata(hFig, 'hSV3DSet');	

	if isempty(hSV3DSet)
		hSV3DSet = 1;
	end

	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		% First, set it to the first "set" of images
% 		hSurf = getappdata(hSV3D, 'hSurf');
% 		setappdata(hSV3D, 'lastCurrentObj', hSurf(1));

		fcnHandles = getappdata(hSV3D, 'fcnHandles');
		fcnHandles.SetTimeSlice(hSV3D, newTime, newSlice, hSV3DSet)
	end

	% Apply to another SimpleViewer, if applicable
	hSV = getappdata(hFig, 'hSV');
	if ~isempty(hSV) && ishandle(hSV)
		fcnHandles = getappdata(hSV, 'fcnHandles');
		fcnHandles.SetTimeSlice(hSV, newTime, newSlice)
	end

	lastCallTime = clock;
% 	pause(1/30)
end

%% -----------------------------------------------------------------------------
function SetTimeSliceImageOnly(hFig, newTime, newSlice)
% Sets the time/slice for all subplots within a figure
% 	if isMultipleCall(); return; end
	prefs       = getappdata(hFig,        'prefs');
	persistent lastCallTime
	if ~isempty(lastCallTime)
		if etime(clock, lastCallTime) < 1/prefs.maxFrameRate
% 			disp('Exeeding frame rate limiter!')
			return
		end
	end

	currentAxes = get(hFig,               'CurrentAxes');
	hAxes       = getappdata(hFig,        'hAxes');
% 	prefs       = getappdata(hFig,        'prefs');
	bManualZoom = getappdata(hFig,        'bManualZoom');
	imgData     = getappdata(currentAxes, 'imgData');
	oldTime     = getappdata(currentAxes, 'currentTime');
	oldSlice    = getappdata(currentAxes, 'currentSlice');
	oldTimeImg  = getappdata(currentAxes, 'currentTimeImg');
	oldSliceImg = getappdata(currentAxes, 'currentSliceImg');

	if isempty(imgData)
		return
	end

	% oldTimeImg and oldSliceImg represent the currently displayed time/slice if
	% different from the current time/slice
	if ~isempty(oldTimeImg)
		oldTime = oldTimeImg;
	end

	if ~isempty(oldSliceImg)
		oldSlice = oldSliceImg;
	end

% 	if (newTime == oldTime) && (newSlice == oldSlice)
% 		return
% 	end

% 	% ROI drawing places some limitations on what you can change
% 	if ~isempty(getappdata(currentAxes, 'inProgressRoi'))
% 		if newSlice ~= oldSlice
% 			set(getappdata(hFig, 'hStatusText'), 'String', 'Error: Cannot change slice while drawing an ROI', 'Color', 'r');
% 			return
% 		end
% 		if ~getappdata(hFig, 'staticRoi') && (newTime ~= oldTime)
% 			set(getappdata(hFig, 'hStatusText'), 'String', 'Error: Cannot change time while drawing an ROI if staticRoi is off', 'Color', 'r');
% 			return
% 		end
% 	end

	% Auto-resize means that the window size will change to keep the same zoom ratio as the previous frame
	if prefs.autoResize
		scale = CalculateCurrentScale(hFig);
		if scale < 0.5
			scale = 0.5;  % Don't allow the window to become too small
		end
	end

	if (newSlice ~= oldSlice)
		% When changing slice, "clip" the time if necessary
		newSlice = mod(newSlice-1, length(imgData))+1;
		newTime  = min([length(imgData{newSlice}) newTime]);
	elseif (newTime ~= oldTime)
 		% Slice is already valid
		newTime  = mod(newTime -1, length(imgData{newSlice}))+1;
	end

	% Update all subplots
	for jAxis = reshape(hAxes(hAxes ~= 0), 1, [])
		newSubTime    = newTime;
		newSubSlice   = newSlice;
		imgData       = getappdata(jAxis, 'imgData');
		txtData       = getappdata(jAxis, 'txtData');
		roi           = getappdata(jAxis, 'roi');
		roiHandles    = getappdata(jAxis, 'roiHandles');
		hImage        = getappdata(jAxis, 'hImage');
		hText         = getappdata(jAxis, 'hText');
		inProgressRoi = getappdata(jAxis, 'inProgressRoi');

		if isempty(imgData)
			continue
		end

		% Range check (these are clipped, not cyclic)
		newSubSlice = min([newSubSlice length(imgData             )]);
		newSubTime  = min([newSubTime  length(imgData{newSubSlice})]);

% 		if (newSubSlice == getappdata(jAxis, 'currentSlice')) && (newSubTime  == getappdata(jAxis, 'currentTime'))
% 			continue
% 		end

		% Image
		oldImgSize = size(get(hImage, 'CData'));
		newImgSize = size(imgData{newSubSlice}{newSubTime});
		set(hImage, 'CData',  imgData{newSubSlice}{newSubTime});

		if isempty(bManualZoom) || ~bManualZoom
			if ~isempty(imgData{newSubSlice}{newSubTime}) && ((oldImgSize(1) ~= newImgSize(1)) || (oldImgSize(2) ~= newImgSize(2)))
				set(jAxis, 'XLim', [0.5 newImgSize(2)+0.5])
				set(jAxis, 'YLim', [0.5 newImgSize(1)+0.5])
			end
		end

		% Text
% 		strText = ComposeText(txtData{newSubSlice}{min([newSubTime  length(txtData{newSubSlice})])}, ...
% 													newSubSlice, length(imgData), ...
% 													newSubTime,  length(imgData{newSubSlice}));
% 		set(hText, 'String', strText);

		% This has to be done before calling UpdateRoiHandles()
		setappdata(jAxis, 'currentSliceImg', newSubSlice);
		setappdata(jAxis, 'currentTimeImg',  newSubTime );

% 		if ~isempty(inProgressRoi)
% 			continue
% 		end
% 
% 		% Update ROI handles for rois on this time/slice
% 		for iRoi = 1:length(roi{newSubSlice}{newSubTime})
% % 			% TODO: this makes it slower,  but need this to update roiPerimeter/roiArea
% % 			UpdateRoiData(hFig, iRoi);
% 			UpdateRoiHandles(jAxis, iRoi)
% 		end
% 
% 		% Hide any additional plot handles
% 		for iRoi = length(roi{newSubSlice}{newSubTime})+1:length(roiHandles)
% 			set(roiHandles(iRoi).line,     'XData', [],   'YData', []);
% 			set(roiHandles(iRoi).vertices, 'XData', [],   'YData', []);
% 			set(roiHandles(iRoi).center,   'XData', [],   'YData', []);
% 			if ~strcmp(prefs.roiLabelMode, 'none')
% 				set(roiHandles(iRoi).label,    'Position', nan(1,2));
% 			end
% 		end
	end

% 	setappdata(currentAxes, 'lastTime',  oldTime);
% 	setappdata(currentAxes, 'lastSlice', oldSlice);

	% New auto window level mode
	bAutoWindowLevel = getappdata(hFig, 'bAutoWindowLevel');
	if ~isempty(bAutoWindowLevel) && bAutoWindowLevel
		AutoWindowLevel(hFig);
	end
	
	% FIXME: This should be a remnant of when drawing points and it required changing the current axis
% 	set(hFig, 'CurrentAxes', currentAxes);
	if prefs.autoResize
		SetFigureScale(hFig, scale);
	end
% 	drawnow expose
% 	pause(1)

% 	if (newSlice ~= oldSlice)
		UpdateRoiFigure(hFig);
% 	end

	% Apply to SimpleViewer3D, if applicable
	hSV3D    = getappdata(hFig, 'hSV3D');
	hSV3DSet = getappdata(hFig, 'hSV3DSet');	

	if isempty(hSV3DSet)
		hSV3DSet = 1;
	end

% % 	if ~isempty(hSV3D) && ishandle(hSV3D)
% % % 		% First, set it to the first "set" of images
% % % 		hSurf = getappdata(hSV3D, 'hSurf');
% % % 		setappdata(hSV3D, 'lastCurrentObj', hSurf(1));
% % 
% % 		fcnHandles = getappdata(hSV3D, 'fcnHandles');
% % 		fcnHandles.SetTimeSliceImageOnly(hSV3D, newTime, newSlice, hSV3DSet)
% % 	end
% % 
% % 	% Apply to another SimpleViewer, if applicable
% % 	hSV = getappdata(hFig, 'hSV');
% % 	if ~isempty(hSV) && ishandle(hSV)
% % 		fcnHandles = getappdata(hSV, 'fcnHandles');
% % 		fcnHandles.SetTimeSliceImageOnly(hSV, newTime, newSlice)
% % 	end

	lastCallTime = clock;
% 	pause(1/30)
end

%%
function MarkFrame(hFig)
	currentAxes  = get(hFig,               'CurrentAxes');
% 	hAxes        = getappdata(hFig,        'hAxes');
	imgData      = getappdata(currentAxes, 'imgData');
	currentTime  = getappdata(currentAxes, 'currentTime');
	currentSlice = getappdata(currentAxes, 'currentSlice');
	markedData   = getappdata(currentAxes, 'markedData');

	if isempty(markedData)
		markedData = false(numel(imgData), numel(imgData{currentSlice}));
	end
	
	markedData(currentSlice,currentTime) = ~markedData(currentSlice,currentTime);
	
	setappdata(currentAxes, 'markedData', markedData);
end
%% -----------------------------------------------------------------------------
% function MakeCircle(hFig, evt)
% % Creates a "force field" circle around the mouse cursor that displaces vertices
% 
% 	% Update the selected subplot, ROI
% 	MouseDownSelect(hFig, evt);
% 
% 	% If the click occurs on an existing ROI object, abort and call that object's callback
% 	hObj  = get(hFig, 'CurrentObject');
% 	hFcn  = get(hObj, 'ButtonDownFcn');
% 	if ~isempty(hFcn)
% 		MouseDownSelect(hFig, evt);
% 		hFcn(hFig, evt);
% 		return
% 	end
% 
% 	hAxis        = get(hFig,        'CurrentAxes');
% 	staticRoi    = getappdata(hFig, 'staticRoi');
% 	prefs        = getappdata(hFig, 'prefs');
% 
% 	roi          = getappdata(hAxis, 'roi');
% 	roiHandles   = getappdata(hAxis, 'roiHandles');
% 	currentSlice = getappdata(hAxis, 'currentSlice');
% 	currentTime  = getappdata(hAxis, 'currentTime');
% 
% 	SelectSubplot(hFig, hAxis);
% % 	switch(get(hFig, 'SelectionType'))
% % 		case 'normal'
% % 	end
% 	
% 	% Get the circle data
% 	hCircle = getappdata(hAxis, 'hCircle');
% 	XData = get(hCircle, 'XData');
% 	YData = get(hCircle, 'YData');	
% 
% 	% Next ROI is the first empty slot (i.e. if ROI 2 of 3 is deleted, the next ROI drawn is index 2)
% 	if isempty(roi{currentSlice}{currentTime})
% 		iRoi = [];
% 	else
% 		iRoi = find(cellfun(@(x) isempty(x), roi{currentSlice}{currentTime}), 1);
% 	end
% 
% 	if isempty(iRoi)
% 		iRoi = length(roi{currentSlice}{currentTime})+1;
% 	end
% 
% 	roi{currentSlice}{currentTime}{iRoi} = [XData; YData]';
% 
% 	% For staticRoi, propagate it to all other time frames
% 	if staticRoi
% 		for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 			roi{currentSlice}{iTime}{iRoi} = roi{currentSlice}{currentTime}{iRoi};
% 		end
% 	end
% 
% 	% Create plot objects if necessary
% 	if isempty(roiHandles)
% 		roiHandles = struct('line', [],  'vertices', [],  'center', [],  'label', []);
% 	end
% 
% 	if iRoi > length(roiHandles)
% 		roiHandles(iRoi) = struct('line', [],  'vertices', [],  'center', [],  'label', []);
% 	end
% 
% 	roiColor = prefs.roiColors(mod(iRoi-1, length(prefs.roiColors))+1);
% 	if isempty(roiHandles(iRoi).line)
% 		roiHandles(iRoi).line   = plot(0, 0, ...
% 																	 'Color',      roiColor, ...
% 																	 'LineStyle',  '-', ...
% 																	 'Marker',     'none');
% 	end
% 
% % 	if isempty(roiHandles(iRoi).vertices)
% 	for iVertex = 1:size(roi{currentSlice}{currentTime}{iRoi},1)
% 		if iVertex > numel(roiHandles(iRoi).vertices)
% 			roiHandles(iRoi).vertices(iVertex) = plot(0, 0, ...
% 			                                          'Color',      roiColor, ...
% 			                                          'LineStyle',  'none', ...
% 			                                          'Marker',     '.', ...
% 			                                          'MarkerSize', 6);
% 		end
% 	end
% 
% 	if isempty(roiHandles(iRoi).center)
% 		roiHandles(iRoi).center = plot(0, 0, ...
% 																	 'Color',      roiColor, ...
% 																	 'LineStyle',  'none', ...
% 																	 'Marker',     '+', ...
% 																	 'MarkerSize', 6);
% 	end
% 
% 	if isempty(roiHandles(iRoi).label) && ~prefs.hideRoiLabels
% 		roiHandles(iRoi).label  = text(0, 0, num2str(iRoi), ...
% 																	 'Color',      roiColor, ...
% 																	 'HitTest',    'off');
% 	end
% 
% 	% Set positions for ROI's plot objects
% 	if ~strcmp(prefs.interpAlgorithm, 'linear')
% 		tmpInterpRoi = CubicInterp(roi{currentSlice}{currentTime}{iRoi}, prefs.interpNPts, prefs.interpAlgorithm);
% 		set(roiHandles(iRoi).line, 'XData', tmpInterpRoi(:,1), 'YData', tmpInterpRoi(:,2));
% 	else
% 		set(roiHandles(iRoi).line, 'XData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],1), ...
% 		                           'YData', roi{currentSlice}{currentTime}{iRoi}([1:end 1],2));
% 	end
% 
% 	for iVertex = 1:size(roi{currentSlice}{currentTime}{iRoi},1)
% 		set(roiHandles(iRoi).vertices(iVertex), 'XData', roi{currentSlice}{currentTime}{iRoi}(iVertex,1), ...
% 		                                        'YData', roi{currentSlice}{currentTime}{iRoi}(iVertex,2));
% 	end
% 
% 	set(roiHandles(iRoi).center,   'XData', mean(roi{currentSlice}{currentTime}{iRoi}(:,1)), ...
% 																 'YData', mean(roi{currentSlice}{currentTime}{iRoi}(:,2)));
% 
% 	if ~prefs.hideRoiLabels
% 		set(roiHandles(iRoi).label,    'Position', [min(roi{currentSlice}{currentTime}{iRoi}(:,1))-3 min(roi{currentSlice}{currentTime}{iRoi}(:,2))-3 0]);
% 	end
% 
% 	set(roiHandles(iRoi).center,   'ButtonDownFcn', @StartRoiMove);
% 	set(roiHandles(iRoi).vertices, 'ButtonDownFcn', @VertexButtonDown);
% 	set(roiHandles(iRoi).line,     'ButtonDownFcn', @AddVertexFromEdge);
% 
% 	% Save back to figure
% 	setappdata(hAxis, 'roi',           roi);
% 	setappdata(hAxis, 'roiHandles',    roiHandles);
% 	UpdateRoiData(hFig, iRoi);
% 
% 	% Apply to SimpleViewer3D, if applicable
% 	hSV3D = getappdata(hFig, 'hSV3D');
% 	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		hAxes = getappdata(hFig, 'hAxes');
% 		if (hAxis == hAxes(1))
% 			setappdata(hSV3D, 'roi', roi);
% 			fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 			fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 		end
% 	end
% end

function MakeCircle(hFig, evt)
% Creates a "force field" circle around the mouse cursor that displaces vertices

	% These shouldn't actually be necessary.  The circle's move function updates
	% the subplot, and we don't want to be able to click on any ROI bits
% 	% Update the selected subplot, ROI
% 	MouseDownSelect(hFig, evt);
% 
% 	% If the click occurs on an existing ROI object, abort and call that object's callback
% 	% TODO: Consider a different approach so that we can draw concentric circles
% 	hObj  = get(hFig, 'CurrentObject');
% 	hFcn  = get(hObj, 'ButtonDownFcn');
% 	if ~isempty(hFcn)
% 		MouseDownSelect(hFig, evt);
% 		hFcn(hFig, evt);
% 		return
% 	end

	CurrentAxes  = get(hFig,               'CurrentAxes');
	roi          = getappdata(CurrentAxes, 'roi');
	currentSlice = getappdata(CurrentAxes, 'currentSlice');
	currentTime  = getappdata(CurrentAxes, 'currentTime');

	SelectSubplot(hFig, CurrentAxes);
% 	switch(get(hFig, 'SelectionType'))
% 		case 'normal'
% 	end
	
	% Get the circle data
	hCircle = getappdata(CurrentAxes, 'hCircle');
	XData   = get(       hCircle,     'XData');
	YData   = get(       hCircle,     'YData');	

	% Take out the last point (it's duplicated)
	XData = XData(1:end-1);
	YData = YData(1:end-1);

	% Next ROI is the first empty slot (i.e. if ROI 2 of 3 is deleted, the next ROI drawn is index 2)
	if isempty(roi{currentSlice}{currentTime})
		iRoi = [];
	else
		iRoi = find(cellfun(@(x) isempty(x), roi{currentSlice}{currentTime}), 1);
	end

	if isempty(iRoi)
		iRoi = length(roi{currentSlice}{currentTime})+1;
	end

	roi{currentSlice}{currentTime}{iRoi} = [XData; YData]';
	setappdata(CurrentAxes, 'roi', roi);
	UpdateRoiData(hFig, iRoi);

	DoStaticAndTranscribeAnd3D(hFig, iRoi)
end

function MoveCircle(hFig, evt)
	CurrentAxes  = get(hFig, 'CurrentAxes');
	CurrentPoint = get(CurrentAxes, 'CurrentPoint');

	hCircle = getappdata(CurrentAxes, 'hCircle');
	prefs = getappdata(hFig, 'prefs');
	aspectRatio = prefs.PixelSpacing(2) / prefs.PixelSpacing(1);

	% Check if the current point is outside of the current axes
	XLim = get(CurrentAxes, 'XLim');
	YLim = get(CurrentAxes, 'YLim');

	if (CurrentPoint(1,1) < XLim(1)) || (CurrentPoint(1,1) > XLim(2)) || ...
	   (CurrentPoint(1,2) < YLim(1)) || (CurrentPoint(1,2) > YLim(2))
		hNewAxes = overobj2('type', 'axes', 'hittest', 'on');
		% TODO: Check to see if the new axis is still a child of the figure

		if (get(hNewAxes, 'Parent') ~= hFig)
			return
		end

		if (numel(hNewAxes) == 1) && ishandle(hNewAxes) && (hNewAxes ~= CurrentAxes)
			% Move circle to new subplot
			set(hCircle, 'Parent', hNewAxes)
			setappdata(hNewAxes, 'hCircle', hCircle);
			setappdata(CurrentAxes, 'hCircle', []);
			set(hFig, 'CurrentAxes', hNewAxes)
			SelectSubplot(hFig, hNewAxes);
			MoveCircle(hFig, evt);
			return
		end
	end

	x = prefs.circleRadius*cos(linspace(0,2*pi,prefs.circlePoints)) / aspectRatio + CurrentPoint(1,1);
	y = prefs.circleRadius*sin(linspace(0,2*pi,prefs.circlePoints))               + CurrentPoint(1,2);

	set(hCircle, 'XData', x);
	set(hCircle, 'YData', y);
end

function ScaleCircle(hFig, scale)
	prefs = getappdata(hFig, 'prefs');
	strMode = getappdata(hFig, 'strMode');
	
	switch strMode
		case 'Circle'
			prefs.circleRadius = prefs.circleRadius * scale;
			setappdata(hFig, 'prefs', prefs);
			MoveCircle(hFig);
		case 'Nudge'
			prefs.forceFieldRadius = prefs.forceFieldRadius * scale;
			setappdata(hFig, 'prefs', prefs);
			MoveNudge(hFig, [], false);
	end
end

function AutoMakeCircle(hFig, evt)
% 	aSearchRadii = [5 15];  % [min max] radius [px]
	scaleFactor  = 0.6;     % Scale detected radius by this factor

	aSearchRadii = (11:10:51);  % Mean radii to search [px]

	hAxis        = get(hFig,         'CurrentAxes');
	currentSlice = getappdata(hAxis, 'currentSlice');
	currentTime  = getappdata(hAxis, 'currentTime');
	imgData      = getappdata(hAxis, 'imgData');
	roi          = getappdata(hAxis, 'roi');
	prefs        = getappdata(hFig,  'prefs');

	% Use window leveled image to assist circle detection
	CLim         = get(hAxis,        'CLim');
	img          = WindowLevel(imgData{currentSlice}{currentTime}, diff(CLim), mean(CLim));

	% Span of radii to search for
	if numel(aSearchRadii) == 1
		diffRadii = 10;
	else
		diffRadii = diff(aSearchRadii(1:2));
	end

	% Search for circles using the Hough transform
	for i = 1:numel(aSearchRadii)
		[cCenters{i}, cRadii{i}] = imfindcircles(interp2(img), aSearchRadii(i) + diffRadii*[-1 1]/2);
	end

	% Find the mean radius that has the most results and use that
	nCircles = cellfun(@(x) numel(x), cRadii);
	indGood = find(nCircles == max(nCircles),1);

	aCenters = cCenters{indGood};
	aRadii   = cRadii{  indGood};

	% Used interp2 for low res images, so shift and scale accordingly
	aCenters = (aCenters-1)./2 + 1;
	aRadii   = (aRadii  -1)./2 + 1;

	% Be a little more conservative with the radius and make them all the same
	aRadii = mean(aRadii) * scaleFactor * ones(size(aRadii));

	% Manually identify the order
	hCircles = viscircles(aCenters, aRadii);
	i = 1;
	while(1)
		[x, y] = ginput(1);
		% "Enter" results in an empty result, so break out of loop
		if isempty(x)
			break
		end
		xClick(i) = x;
		yClick(i) = y;
		hClick(i) = text(xClick(i), yClick(i), num2str(i), 'HorizontalAlignment', 'center', 'color', 'r');
		i = i+1;
	end
	delete(hCircles)
	if exist('hClick', 'var')
		delete(hClick)
	end
	set(hFig, 'currentch', char(1));  % reset 'CurrentCharacter'

	if ~exist('xClick', 'var') || isempty(xClick)
		return
	end

	% Find the Hough circle closest to each clicked point and add to ROI list
	aNewRois = [];
	for i = 1:numel(xClick)
		dist     = sum((aCenters - repmat([xClick(i) yClick(i)], [size(aCenters,1) 1])).^2, 2);
		indFound = find(dist == min(dist), 1);

		x = aRadii(indFound)*cos(linspace(0,2*pi,prefs.circlePoints)) + aCenters(indFound,1);
		y = aRadii(indFound)*sin(linspace(0,2*pi,prefs.circlePoints)) + aCenters(indFound,2);

		tmpRoi = [x; y]';

		% Next ROI is the first empty slot (i.e. if ROI 2 of 3 is deleted, the next ROI drawn is index 2)
		if isempty(roi{currentSlice}{currentTime})
			iRoi = [];
		else
			iRoi = find(cellfun(@(x) isempty(x), roi{currentSlice}{currentTime}), 1);
		end

		if isempty(iRoi)
			iRoi = length(roi{currentSlice}{currentTime})+1;
		end

		roi{currentSlice}{currentTime}{iRoi} = tmpRoi;
		aNewRois = cat(2, aNewRois, iRoi);
	end
	setappdata(hAxis, 'roi', roi);

	% UpdateRoiData really needs to be re-written to support multiple ROIs...
	for jRoi = aNewRois
		UpdateRoiData(hFig, jRoi);
	end
	DoStaticAndTranscribeAnd3D(hFig, aNewRois)
end

%% -----------------------------------------------------------------------------
function StartNudge(hFig, evt)
% Creates a "force field" circle around the mouse cursor that displaces vertices

	% Update the selected subplot, ROI
	hAxis = get(hFig, 'CurrentAxes');
	SelectSubplot(hFig, hAxis);

	% Abort if already nudging (i.e. right click while left click is still held)
	if ~isempty(getappdata(hAxis, 'hForceField'))
		return
	end

% Note: Disabling this here only means that the other object's callback isn't
% called before the rest of this function.  It still runs AFTER!

% 	% If the click occurs on an existing ROI object, abort and call that object's callback
% 	hObj  = get(hFig, 'CurrentObject');
% 	hFcn  = get(hObj, 'ButtonDownFcn');
% 	if ~isempty(hFcn)
% 		MouseDownSelect(hFig, evt);
% 		hFcn(hFig, evt);
% 		return
% 	end

	selectedRoi  = getappdata(hFig,  'selectedRoi');
	roi          = getappdata(hAxis, 'roi');
	currentSlice = getappdata(hAxis, 'currentSlice');
	currentTime  = getappdata(hAxis, 'currentTime');
	prefs        = getappdata(hFig,  'prefs');
	CurrentPoint = get(       hAxis, 'CurrentPoint');

	aspectRatio = prefs.PixelSpacing(2) / prefs.PixelSpacing(1);

	% Find the nearest vertex
	if isempty(roi{currentSlice}{currentTime})
		prefs.forceFieldRadius = 25;
	else
		if isempty(selectedRoi)
			selectedRoi = 1:numel(roi{currentSlice}{currentTime});
		end

		aMinDist = realmax*ones(1,max(selectedRoi));
		for jRoi = selectedRoi
			if jRoi > numel(roi{currentSlice}{currentTime})
				continue
			end
			if ~isempty(roi{currentSlice}{currentTime}{jRoi})
				diffDist = roi{currentSlice}{currentTime}{jRoi} - repmat(CurrentPoint(1,1:2), [size(roi{currentSlice}{currentTime}{jRoi},1) 1]);
				diffDist(:,1) = diffDist(:,1) * aspectRatio;
				dist = sqrt(sum((diffDist).^2, 2));
				aMinDist(jRoi) = min(sqrt(sum(diffDist.^2, 2)));
			end
		end
		prefs.forceFieldRadius = min(aMinDist) * 0.85;
	end
	setappdata(hFig, 'prefs', prefs);

	hAxis = get(hFig, 'CurrentAxes');
	switch(get(hFig, 'SelectionType'))
		case 'normal'
% 			set(hFig, 'Pointer', 'crosshair');
			% Draw force field
			CurrentPoint = get(hAxis, 'CurrentPoint');
			x = prefs.forceFieldRadius*cos(0:pi/50:2*pi) / aspectRatio + CurrentPoint(1,1);
			y = prefs.forceFieldRadius*sin(0:pi/50:2*pi)               + CurrentPoint(1,2);

			hForceField = plot(x,y, 'r-');
			setappdata(hAxis, 'hForceField', hForceField);

% 			setappdata(hFig, 'InitPoint', get(hAxis, 'CurrentPoint'))
% 			setappdata(hFig, 'InitCLim',  get(hAxis, 'CLim'))
% 			set(hFig, 'WindowButtonMotionFcn', @MoveNudge);
			set(hFig, 'WindowButtonMotionFcn', @(hFig, evt) MoveNudge(hFig, evt, false));
			set(hFig, 'WindowButtonUpFcn',     @EndNudge);
		case 'alt'
			if isempty(roi{currentSlice}{currentTime})
				prefs.forceFieldRadius = 25;
			else
				if isempty(selectedRoi)
					selectedRoi = 1:numel(roi{currentSlice}{currentTime});
				end

				aMaxDist = zeros(1,max(selectedRoi));
				for jRoi = selectedRoi
					if ~isempty(roi{currentSlice}{currentTime}{jRoi})
						diffDist = roi{currentSlice}{currentTime}{jRoi} - repmat(CurrentPoint(1,1:2), [size(roi{currentSlice}{currentTime}{jRoi},1) 1]);
						diffDist(:,1) = diffDist(:,1) * aspectRatio;
						dist = sqrt(sum((diffDist).^2, 2));
						aMaxDist(jRoi) = max(sqrt(sum(diffDist.^2, 2)));
					end
				end
				prefs.forceFieldRadius = max(aMaxDist) * 1.15;
			end
			setappdata(hFig, 'prefs', prefs);

			% Draw force field
			CurrentPoint = get(hAxis, 'CurrentPoint');
			x = prefs.forceFieldRadius*cos(0:pi/50:2*pi) / aspectRatio + CurrentPoint(1,1);
			y = prefs.forceFieldRadius*sin(0:pi/50:2*pi)               + CurrentPoint(1,2);

			hForceField = plot(x,y, 'b-');
			setappdata(hAxis, 'hForceField', hForceField);

% 			set(hFig, 'WindowButtonMotionFcn', @MoveNudgeOpposite);
			set(hFig, 'WindowButtonMotionFcn', @(hFig, evt) MoveNudge(hFig, evt, true));
			set(hFig, 'WindowButtonUpFcn',     @EndNudge);
	end
	
	% Save the current ROI setfor undo
	setappdata(hAxis, 'prevRoi', roi);
end

function EndNudge(hFig, evt)
% Disperses nudge force field
% 	set(hFig, 'Pointer', 'arrow');

	% Delete displayed circle
	CurrentAxes = get(hFig, 'CurrentAxes');
	hForceField = getappdata(CurrentAxes, 'hForceField');
	delete(hForceField);
	setappdata(CurrentAxes, 'hForceField', []);

	% Automatically interpolate vertices
	prefs = getappdata(hFig, 'prefs');
	if prefs.autoInterpolateRois || strcmp(prefs.interpRois, 'always')
		InterpolateRoi(CurrentAxes);
	end
	
	set(hFig, 'WindowButtonMotionFcn', '');
	set(hFig, 'WindowButtonUpFcn',     '');
end

% function MoveNudge(hFig, evt)
% 	CurrentAxes  = get(hFig, 'CurrentAxes');
% 	currentSlice = getappdata(CurrentAxes, 'currentSlice');
% 	currentTime  = getappdata(CurrentAxes, 'currentTime');
% 	staticRoi    = getappdata(hFig, 'staticRoi');
% 
% 	hForceField = getappdata(CurrentAxes, 'hForceField');
% 	prefs = getappdata(hFig, 'prefs');
% 
% 	aspectRatio = prefs.PixelSpacing(2) / prefs.PixelSpacing(1);
% 	
% 	% Update displayed circle
% 	CurrentPoint = get(CurrentAxes, 'CurrentPoint');
% 	x = prefs.forceFieldRadius*cos(0:pi/25:2*pi) / aspectRatio + CurrentPoint(1,1);
% 	y = prefs.forceFieldRadius*sin(0:pi/25:2*pi)               + CurrentPoint(1,2);
% 
% 	set(hForceField, 'XData', x);
% 	set(hForceField, 'YData', y);
% 	
% 	% Check for vertices that lie within the force field radius
% 	roi = getappdata(CurrentAxes, 'roi');
% 
% 	selectedRoi  = getappdata(hFig,  'selectedRoi');
% 	if isempty(selectedRoi)
% 		selectedRoi = 1:numel(roi{currentSlice}{currentTime});
% 	end
% 
% 	for iRoi = selectedRoi
% 		if isempty(roi{currentSlice}{currentTime}{iRoi})
% 			continue
% 		end
% 
% 		diffDist = roi{currentSlice}{currentTime}{iRoi} - repmat(CurrentPoint(1,1:2), [size(roi{currentSlice}{currentTime}{iRoi},1) 1]);
% 		diffDist(:,1) = diffDist(:,1) * aspectRatio;
% 		dist = sqrt(sum((diffDist).^2, 2));
% 
% 		% TODO: this should be able to run without the for loop
% 		badVertices = find(dist < prefs.forceFieldRadius)';
% 		for iVertex = badVertices
% 			roi{currentSlice}{currentTime}{iRoi}(iVertex,:) = (roi{currentSlice}{currentTime}{iRoi}(iVertex,:) - CurrentPoint(1,1:2)) * prefs.forceFieldRadius/dist(iVertex) + CurrentPoint(1,1:2);
% 
% 			% For staticRoi, propagate it to all other time frames
% 			if staticRoi
% 				for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 					roi{currentSlice}{iTime}(iRoi) = roi{currentSlice}{currentTime}(iRoi);
% 				end
% 			end
% 
% 			% TODO: Can we move this outside the for loop?
% 			setappdata(CurrentAxes, 'roi', roi);
% 
% 			% TODO: Need UpdateRoiData to update area/perimeter, but its current implementation is broken
% % 			UpdateRoiData(hFig, iRoi);
% 			UpdateRoiHandles(CurrentAxes, iRoi);
% 		end
% 	end
% 	
% 	% Apply to SimpleViewer3D, if applicable
% 	hSV3D = getappdata(hFig, 'hSV3D');
% 	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		hAxes = getappdata(hFig, 'hAxes');
% 		if (CurrentAxes == hAxes(1))
% 			setappdata(hSV3D, 'roi', roi);
% 			fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 			fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 		end
% 	end
% end

% function MoveNudgeOpposite(hFig, evt)
% 	CurrentAxes  = get(hFig, 'CurrentAxes');
% 	currentSlice = getappdata(CurrentAxes, 'currentSlice');
% 	currentTime  = getappdata(CurrentAxes, 'currentTime');
% 	staticRoi    = getappdata(hFig, 'staticRoi');
% 
% 	hForceField = getappdata(CurrentAxes, 'hForceField');
% 	prefs = getappdata(hFig, 'prefs');
% 
% 	aspectRatio = prefs.PixelSpacing(2) / prefs.PixelSpacing(1);
% 	
% 	% Update displayed circle
% 	CurrentPoint = get(CurrentAxes, 'CurrentPoint');
% 	x = prefs.forceFieldRadius*cos(0:pi/25:2*pi) / aspectRatio + CurrentPoint(1,1);
% 	y = prefs.forceFieldRadius*sin(0:pi/25:2*pi)               + CurrentPoint(1,2);
% 
% 	set(hForceField, 'XData', x);
% 	set(hForceField, 'YData', y);
% 	
% 	% Check for vertices that lie within the force field radius
% 	roi = getappdata(CurrentAxes, 'roi');
% 
% 	selectedRoi  = getappdata(hFig,  'selectedRoi');
% 	if isempty(selectedRoi)
% 		selectedRoi = 1:numel(roi{currentSlice}{currentTime});
% 	end
% 
% 	for iRoi = selectedRoi
% 		if isempty(roi{currentSlice}{currentTime}{iRoi})
% 			continue
% 		end
% 
% 		diffDist = roi{currentSlice}{currentTime}{iRoi} - repmat(CurrentPoint(1,1:2), [size(roi{currentSlice}{currentTime}{iRoi},1) 1]);
% 		diffDist(:,1) = diffDist(:,1) * aspectRatio;
% 		dist = sqrt(sum((diffDist).^2, 2));
% 
% 		% TODO: this should be able to run without the for loop
% 		badVertices = find(dist > prefs.forceFieldRadius)';
% 		for iVertex = badVertices
% 			roi{currentSlice}{currentTime}{iRoi}(iVertex,:) = (roi{currentSlice}{currentTime}{iRoi}(iVertex,:) - CurrentPoint(1,1:2)) * prefs.forceFieldRadius/dist(iVertex) + CurrentPoint(1,1:2);
% 			
% 			% For staticRoi, propagate it to all other time frames
% 			if staticRoi
% 				for iTime = [1:currentTime-1 currentTime+1:length(roi{currentSlice})]
% 					roi{currentSlice}{iTime}(iRoi) = roi{currentSlice}{currentTime}(iRoi);
% 				end
% 			end
% 
% 			setappdata(CurrentAxes, 'roi', roi);
% % 			UpdateRoiData(hFig, iRoi);
% 			UpdateRoiHandles(CurrentAxes, iRoi);
% 		end
% 	end
% 	
% 	% Apply to SimpleViewer3D, if applicable
% 	hSV3D = getappdata(hFig, 'hSV3D');
% 	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		hAxes = getappdata(hFig, 'hAxes');
% 		if (CurrentAxes == hAxes(1))
% 			setappdata(hSV3D, 'roi', roi);
% 			fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 			fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 		end
% 	end
% end

function MoveNudge(hFig, evt, bOpposite)
	CurrentAxes  = get(hFig,               'CurrentAxes');
	currentSlice = getappdata(CurrentAxes, 'currentSlice');
	currentTime  = getappdata(CurrentAxes, 'currentTime');
	hForceField  = getappdata(CurrentAxes, 'hForceField');
	prefs = getappdata(hFig, 'prefs');

	aspectRatio = prefs.PixelSpacing(2) / prefs.PixelSpacing(1);
	
	% Update displayed circle
	CurrentPoint = get(CurrentAxes, 'CurrentPoint');
	x = prefs.forceFieldRadius*cos(0:pi/25:2*pi) / aspectRatio + CurrentPoint(1,1);
	y = prefs.forceFieldRadius*sin(0:pi/25:2*pi)               + CurrentPoint(1,2);

	set(hForceField, 'XData', x);
	set(hForceField, 'YData', y);
	
	% Check for vertices that lie within the force field radius
	roi = getappdata(CurrentAxes, 'roi');

	selectedRoi  = getappdata(hFig,  'selectedRoi');
	if isempty(selectedRoi)
		selectedRoi = 1:numel(roi{currentSlice}{currentTime});
	end

	aChangedRois = zeros(size(selectedRoi));
	for iRoi = 1:numel(selectedRoi)
		jRoi = selectedRoi(iRoi);
		if isempty(roi{currentSlice}{currentTime}{jRoi})
			continue
		end

		diffDist = roi{currentSlice}{currentTime}{jRoi} - repmat(CurrentPoint(1,1:2), [size(roi{currentSlice}{currentTime}{jRoi},1) 1]);
		diffDist(:,1) = diffDist(:,1) * aspectRatio;
		dist = sqrt(sum((diffDist).^2, 2));

		if bOpposite
			badVertices = find(dist > prefs.forceFieldRadius)';
		else
			badVertices = find(dist < prefs.forceFieldRadius)';
		end

		% TODO: this should be able to run without the for loop?
		for iVertex = badVertices
			roi{currentSlice}{currentTime}{jRoi}(iVertex,:) = (roi{currentSlice}{currentTime}{jRoi}(iVertex,:) - CurrentPoint(1,1:2)) * prefs.forceFieldRadius/dist(iVertex) + CurrentPoint(1,1:2);
		end

		if (numel(badVertices) > 0)
			aChangedRois(iRoi) = 1;
		end
	end
	setappdata(CurrentAxes, 'roi', roi);

	% TODO: Need UpdateRoiData to update area/perimeter, but its current implementation is broken
	DoStaticAndTranscribeAnd3D(hFig, selectedRoi(aChangedRois ~= 0));
end

function UndoRoi(hFig)
	CurrentAxes  = get(hFig,               'CurrentAxes');
	prevRoi      = getappdata(CurrentAxes, 'prevRoi');
	roi          = getappdata(CurrentAxes, 'roi');
	currentSlice = getappdata(CurrentAxes, 'currentSlice');
	currentTime  = getappdata(CurrentAxes, 'currentTime');

	if ~isempty(prevRoi)
		setappdata(CurrentAxes, 'roi', prevRoi);
		DoStaticAndTranscribeAnd3D(hFig, 1:numel(prevRoi{currentSlice}{currentTime}));
	end
% 	
% % 		UpdateRoiData(hFig, iRoi);
% 
% 		for iRoi = 1:numel(prevRoi{currentSlice}{currentTime})
% 			UpdateRoiHandles(CurrentAxes, iRoi);
% 		end
% 	end
% 	
% 	% Apply to SimpleViewer3D, if applicable
% 	hSV3D = getappdata(hFig, 'hSV3D');
% 	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		hAxes = getappdata(hFig, 'hAxes');
% 		if (CurrentAxes == hAxes(1))
% 			setappdata(hSV3D, 'roi', roi);
% 			fcnHandles = getappdata(hSV3D, 'fcnHandles');
% 			fcnHandles.UpdateRoiHandles(hSV3D, iRoi);
% 		end
% 	end
end

%% -----------------------------------------------------------------------------
function StartWindowLevel(hFig, evt)
% Called when a mouse-click is started; dispatches to the appropriate window/level function

	% Update the selected subplot, ROI
	MouseDownSelect(hFig, evt);

	% If the click occurs on an existing ROI object, abort and call that object's callback
	hObj  = get(hFig, 'CurrentObject');
	hFcn  = get(hObj, 'ButtonDownFcn');
	if ~isempty(hFcn)
		MouseDownSelect(hFig, evt);
		hFcn(hFig, evt);
		return
	end

	setappdata(hFig, 'bAutoWindowLevel', false)
	hAxis = get(hFig, 'CurrentAxes');
	switch(get(hFig, 'SelectionType'))
		case 'normal'
			% Left Click: Manual adjustment
			setappdata(hFig, 'InitPoint', get(hAxis, 'CurrentPoint'))
			setappdata(hFig, 'InitCLim',  get(hAxis, 'CLim'))
			set(hFig, 'WindowButtonMotionFcn', @AdjustWindowLevel);
			set(hFig, 'WindowButtonUpFcn',     @EndWindowLevel);
		case 'alt'
			% Right Click: Auto window level for current subplot
			AutoWindowLevel(hFig);
			set(getappdata(hFig, 'hStatusText'), 'String', 'Auto Window/Level');
		case 'open'
			% Double Click: Reset to auto window level for each subplot
			ResetWindowLevel(hFig);
			set(getappdata(hFig, 'hStatusText'), 'String', 'Reset Window/Level');
		case 'extend'
			% Shift Click: Auto window level for each subplot, using only visible data
			AutoWindowLevel(hFig);
			setappdata(hFig, 'bAutoWindowLevel', true)
			set(getappdata(hFig, 'hStatusText'), 'String', 'Dynamic Auto Window/Level');
	end
end

% ------------------------------------------------------------------------------
function EndWindowLevel(hFig, evt)
% Called when a mouse-click is lifted; stops mouse movements from manually adjusting window level

	set(hFig, 'WindowButtonMotionFcn', '');
	set(hFig, 'WindowButtonUpFcn',     '');
end

% ------------------------------------------------------------
function AdjustWindowLevel(hFig, evt)
% Left/right adjusts window, up/down adjusts level
% Adapted from 'imagescn'.
	hCurrentAxis = get(hFig, 'CurrentAxes');

	InitPoint = getappdata(hFig, 'InitPoint');
	InitCLim  = getappdata(hFig, 'InitCLim');
	prefs = getappdata(hFig, 'prefs');

	CurrentPoint = get(hCurrentAxis, 'CurrentPoint');

	XLim = get(hCurrentAxis, 'XLim');
	YLim = get(hCurrentAxis, 'YLim');

	window = (InitCLim(2) - InitCLim(1));
	level  = (InitCLim(2) + InitCLim(1))/2;

	% Use fraction, i.e. relative to position to the originally clicked point
	% to determine the change in window and level
	delta = CurrentPoint(1,1:2) - InitPoint(1,1:2);
	delta = delta .* 5;

	% To change WL sensitivity to position, change exponent to bigger/ smaller odd number 
	sensitivityFactor = 3;
	newLevel =   level  + level  * (delta(2) / diff(YLim))^sensitivityFactor;
	newWindow =  window + window * (delta(1) / diff(XLim))^sensitivityFactor;

	if prefs.disableLevelAdjust
		newLevel = level;
	end
	
	% make sure clims stay ascending
	if (newWindow <= 0)
		newWindow = 0.1;
	end

	aRange = [newLevel - newWindow/2 , newLevel + newWindow/2];
% 	% Set for all axes
% 	hAxes = getappdata(hFig, 'hAxes');
% 	set(hAxes(hAxes ~= 0), 'CLim', [newLevel - newWindow/2 , newLevel + newWindow/2]);
% 
% 	% Apply to m-mode image
% 	hAxesMMode =  getappdata(hFig, 'hAxesMMode');
% 	if ~isempty(hAxesMMode)
% 		for iAxes = 1:numel(hAxesMMode)
% 			set(hAxesMMode(iAxes), 'CLim', [newLevel - newWindow/2 , newLevel + newWindow/2]);
% 		end
% 	end
% 
% 	% Apply to SimpleViewer3D, if applicable
% 	hSV3D = getappdata(hFig, 'hSV3D');
% 	if ~isempty(hSV3D) && ishandle(hSV3D)
% 		set(get(hSV3D, 'CurrentAxes'), 'CLim', [newLevel - newWindow/2 , newLevel + newWindow/2]);
% 	end
		
	hAxes = getappdata(hFig, 'hAxes');
	indCurrSubplot = find(hAxes == hCurrentAxis);

	% Use alt-left-click to window/level only current subplot
	CurrentModifier = get(hFig, 'CurrentModifier');
	if (~isempty(CurrentModifier) && strcmp(CurrentModifier, 'alt'))
		indsGoodAxes = indCurrSubplot;
	elseif iscell(prefs.linkWindowLevel)
		% linkWindowLevel tells which subplots have their window level linked to this
		goodLinkInds = cellfun(@(x) any(x == indCurrSubplot), prefs.linkWindowLevel);
		if any(goodLinkInds)
			indsGoodAxes = prefs.linkWindowLevel{goodLinkInds};
		else
			indsGoodAxes = indCurrSubplot;
		end
	elseif prefs.linkWindowLevel
		indsGoodAxes = find(hAxes ~= 0);
	elseif ~prefs.linkWindowLevel
		indsGoodAxes = indCurrSubplot;
	end

	% Set for all axes
	set(hAxes(indsGoodAxes), 'CLim', aRange);

	% Apply to m-mode image
	hAxesMMode =  getappdata(hFig, 'hAxesMMode');
	if ~isempty(hAxesMMode)
		if ishandle(hAxesMMode(indsGoodAxes))
			set(hAxesMMode(indsGoodAxes), 'CLim', aRange)
		end
	end
	
	% Apply to SimpleViewer3D, if applicable
	hSV3D = getappdata(hFig, 'hSV3D');
	if ~isempty(hSV3D) && ishandle(hSV3D)
		set(get(hSV3D, 'CurrentAxes'), 'CLim', aRange);
	end
end

% ------------------------------------------------------------------------------
function ResetWindowLevel(hFig, evt)
% Each subplot is independently automatically window leveled
	CurrentAxes = get(hFig, 'CurrentAxes');

	% Don't allow this for "centered" images
	prefs = getappdata(hFig, 'prefs');
	if prefs.disableLevelAdjust
		return
	end
	
	hAxes = getappdata(hFig, 'hAxes');
	indCurrSubplot = find(hAxes == CurrentAxes);

	CurrentModifier = get(hFig, 'CurrentModifier');
	if (~isempty(CurrentModifier) && strcmp(CurrentModifier, 'alt'))
		% Use alt-left-click to window/level only current subplot
		indsGoodAxes = indCurrSubplot;
	elseif iscell(prefs.linkWindowLevel)
		% linkWindowLevel tells which subplots have their window level linked to this
		goodLinkInds = cellfun(@(x) any(x == indCurrSubplot), prefs.linkWindowLevel);
		if any(goodLinkInds)
			indsGoodAxes = prefs.linkWindowLevel{goodLinkInds};
		else
			indsGoodAxes = indCurrSubplot;
		end
	elseif prefs.linkWindowLevel
		indsGoodAxes = find(hAxes ~= 0);
	elseif ~prefs.linkWindowLevel
		indsGoodAxes = indCurrSubplot;
	end

	% Set for all axes
	set(hAxes(indsGoodAxes), 'CLimMode', 'auto');

	% Apply to m-mode image
	hAxesMMode =  getappdata(hFig, 'hAxesMMode');
	if ~isempty(hAxesMMode)
		if ishandle(hAxesMMode(indsGoodAxes))
			set(hAxesMMode(indsGoodAxes), 'CLimMode', 'auto')
		end
	end
	
	% Apply to SimpleViewer3D, if applicable
	hSV3D = getappdata(hFig, 'hSV3D');
	if ~isempty(hSV3D) && ishandle(hSV3D)
		set(get(hSV3D, 'CurrentAxes'), 'CLimMode', 'auto');
	end
end

% ------------------------------------------------------------------------------
function AutoWindowLevel(hFig)
% Calculate best window level for current subplot's data, set for all axes
	CurrentAxes = get(hFig, 'CurrentAxes');
	hImage = getappdata(CurrentAxes, 'hImage');
	CData  = get(hImage, 'CData');
	prefs = getappdata(hFig, 'prefs');

	if isempty(CData)
		return
	end

	% Also can't window level an TrueColor image
	if (size(CData, 3) == 3)
		return
	end
	
	XLim = get(CurrentAxes, 'XLim');  XLim(1) = max([1 ceil(XLim(1))]);  XLim(2) = min([size(CData,2) floor(XLim(2))]);
	YLim = get(CurrentAxes, 'YLim');  YLim(1) = max([1 ceil(YLim(1))]);  YLim(2) = min([size(CData,2) floor(YLim(2))]);
	
	% This may fail if we change slice/time and the image has a different size,
	% but XLim/YLim are set for the previous image
	try
		aRange = double(AutoLims(CData(YLim(1):YLim(2), XLim(1):XLim(2))));
	catch
		aRange = double(AutoLims(CData));
	end

	% Abort if the new range is zero
	if all(aRange == 0) || all(isnan(aRange)) || (diff(aRange) == 0)
		return
	end

	% Reset the original level for "centered" images
	if prefs.disableLevelAdjust
		CLim  = get(CurrentAxes, 'CLim');
		level = (CLim(2) + CLim(1))/2;

		halfWindow = max(abs(aRange - level));
		aRange = level + halfWindow*[-1 1];
	end

	hAxes = getappdata(hFig, 'hAxes');
	indCurrSubplot = find(hAxes == CurrentAxes);

	CurrentModifier = get(hFig, 'CurrentModifier');
	if (~isempty(CurrentModifier) && strcmp(CurrentModifier, 'alt'))
		% Use alt-left-click to window/level only current subplot
		indsGoodAxes = indCurrSubplot;
	elseif iscell(prefs.linkWindowLevel)
		% linkWindowLevel tells which subplots have their window level linked to this
		goodLinkInds = cellfun(@(x) any(x == indCurrSubplot), prefs.linkWindowLevel);
		if any(goodLinkInds)
			indsGoodAxes = prefs.linkWindowLevel{goodLinkInds};
		else
			indsGoodAxes = indCurrSubplot;
		end
	elseif prefs.linkWindowLevel
		indsGoodAxes = find(hAxes ~= 0);
	elseif ~prefs.linkWindowLevel
		indsGoodAxes = indCurrSubplot;
	end

	% Set for all axes
	set(hAxes(indsGoodAxes), 'CLim', aRange);

	% Apply to m-mode image
	hAxesMMode =  getappdata(hFig, 'hAxesMMode');
	if ~isempty(hAxesMMode)
		if ishandle(hAxesMMode(indsGoodAxes))
			set(hAxesMMode(indsGoodAxes), 'CLim', aRange)
		end
% 		for iAxes = 1:numel(hAxesMMode)
% 			set(hAxesMMode(iAxes), 'CLim', aRange)
% 		end
	end
	
	% Apply to SimpleViewer3D, if applicable
	hSV3D = getappdata(hFig, 'hSV3D');
	if ~isempty(hSV3D) && ishandle(hSV3D)
		set(get(hSV3D, 'CurrentAxes'), 'CLim', aRange);
	end
end

%% -----------------------------------------------------------------------------
function ToggleFreezeColormap(hAxis)
% - (Un)freeze the colormap and window level of a subplot

	imgDataIndexed = getappdata(hAxis, 'imgDataIndexed');
	
	if isempty(imgDataIndexed)
		FreezeColormap(hAxis);
	else
		UnFreezeColormap(hAxis);
	end
end

% -----------------------------------------------------------------------------
function FreezeColormap(hAxis)
% - Freeze the colormap and window level of a subplot

	hFig           = get(       hAxis, 'Parent');
	CLim           = get(       hAxis, 'CLim');
	CMap           = get(       hFig,  'Colormap');
	imgData        = getappdata(hAxis, 'imgData');
	imgDataIndexed = getappdata(hAxis, 'imgDataIndexed');

	% If previously frozen, abort
	if ~isempty(imgDataIndexed)
		return
	end

	% Convert everything to TrueColor indexed
	imgDataTrueColor = cell(size(imgData));
	for iSli = 1:numel(imgData)
		imgDataTrueColor{iSli} = cell(size(imgData{iSli}));
		for iTime = 1:numel(imgData{iSli})
			imgDataTrueColor{iSli}{iTime} = ind2rgb(round(WindowLevel(imgData{iSli}{iTime}, diff(CLim), mean(CLim))*size(CMap,1)), CMap);
		end

		% Make NaNs black (ind2rgb converts them to the last index in CMap)
		nanMap = repmat(~isnan(imgData{iSli}{iTime}), [1 1 3]);
		imgDataTrueColor{iSli}{iTime} = imgDataTrueColor{iSli}{iTime} .* nanMap;
	end

	setappdata(hAxis, 'imgDataIndexed', imgData);
	setappdata(hAxis, 'imgData',        imgDataTrueColor);

	% Update the currently displayed image
	currentSlice = getappdata(hAxis, 'currentSlice');
	currentTime  = getappdata(hAxis, 'currentTime');
	hImage       = getappdata(hAxis, 'hImage');

	set(hImage, 'CData', imgDataTrueColor{currentSlice}{currentTime});
end

% ------------------------------------------------------------------------------
function UnFreezeColormap(hAxis)
% - Unfreeze the colormap and window level of a subplot

% 	hFig           = get(       hAxis, 'Parent');
% 	CLim           = get(       hAxis, 'CLim');
% 	CMap           = get(       hFig,  'Colormap');
% 	imgData        = getappdata(hAxis, 'imgData');
	imgDataIndexed = getappdata(hAxis, 'imgDataIndexed');

	% Unfrozen, or not enough data to un-freeze
	if isempty(imgDataIndexed)
		return
	end

	setappdata(hAxis, 'imgDataIndexed', []);
	setappdata(hAxis, 'imgData',        imgDataIndexed);

	% Update the currently displayed image
	currentSlice = getappdata(hAxis, 'currentSlice');
	currentTime  = getappdata(hAxis, 'currentTime');
	hImage       = getappdata(hAxis, 'hImage');

	set(hImage, 'CData', imgDataIndexed{currentSlice}{currentTime});
end

%% -----------------------------------------------------------------------------
function ToggleCircle(hFig)
	hAxis = get(hFig, 'CurrentAxes');

	if ~isempty(get(hFig, 'WindowButtonDownFcn')) && strcmp(func2str(get(hFig, 'WindowButtonDownFcn')), func2str(@MakeCircle))
		SetNormalMode(hFig);
		return
	end

	setptr(hFig, 'datacursor');
	
	% Draw a temp circle
	prefs = getappdata(hFig, 'prefs');
	aspectRatio = prefs.PixelSpacing(2) / prefs.PixelSpacing(1);

	% Simulate a mouse click to update 'CurrentPoint'.  Unfortunately, it also
	% queues up this mouse click which results in laying down an actual ROI
% 	robot = java.awt.Robot;
% 	robot.mousePress(java.awt.event.InputEvent.BUTTON1_MASK);
% 	robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);
	CurrentPoint = get(hAxis, 'CurrentPoint');

	x = prefs.circleRadius*cos(linspace(0,2*pi,prefs.circlePoints)) / aspectRatio + CurrentPoint(1,1);
	y = prefs.circleRadius*sin(linspace(0,2*pi,prefs.circlePoints))               + CurrentPoint(1,2);

	hCircle = getappdata(hAxis, 'hCircle');
	if ~isempty(hCircle)
		% Recycle an already existing hCircle (save the environment!)
		set(hCircle, 'XData', x)
		set(hCircle, 'YData', y)
	else
		hCircle = plot(x,y, 'r.-', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
		setappdata(hAxis, 'hCircle', hCircle);
	end
	setappdata(hFig, 'strMode', 'Circle');
	set(hFig, 'WindowButtonMotionFcn', @MoveCircle);

	set(hFig, 'WindowButtonDownFcn', @MakeCircle);

	UpdateStatusText(hFig);
end

function ToggleNudgeMode(hFig)
	if ~isempty(get(hFig, 'WindowButtonDownFcn')) && strcmp(func2str(get(hFig, 'WindowButtonDownFcn')), func2str(@StartNudge))
		SetNormalMode(hFig);
		return
	end

	SetNormalMode(hFig);  % Clean up any cruft from other modes (e.g. display from circle tool)
	set(hFig, 'WindowButtonDownFcn',   @StartNudge);

	pData = getappdata(hFig, 'pData');
	set(hFig, 'Pointer', 'custom', 'PointerShapeCData', pData.forceField, 'PointerShapeHotSpot', [8 8])

	setappdata(hFig, 'strMode', 'Nudge');
	UpdateStatusText(hFig);
end

function TogglePanZoomMode(hFig)
	if ~isempty(get(hFig, 'WindowButtonDownFcn')) && strcmp(func2str(get(hFig, 'WindowButtonDownFcn')), func2str(@StartPanZoom))
		SetNormalMode(hFig);
		return
	end

	SetNormalMode(hFig);  % Clean up any cruft from other modes (e.g. display from circle tool)
	set(hFig, 'WindowButtonDownFcn', @StartPanZoom);

	pData = getappdata(hFig, 'pData');
	set(hFig, 'Pointer', 'custom', 'PointerShapeCData', pData.panZoom, 'PointerShapeHotSpot', [8 8])

	setappdata(hFig, 'strMode', 'Pan/Zoom');
	UpdateStatusText(hFig);
end


function ToggleWindowLevelMode(hFig)
	if ~isempty(get(hFig, 'WindowButtonDownFcn')) && strcmp(func2str(get(hFig, 'WindowButtonDownFcn')), func2str(@StartWindowLevel))
		SetNormalMode(hFig);
		return
	end

	SetNormalMode(hFig);  % Clean up any cruft from other modes (e.g. display from circle tool)
	set(hFig, 'WindowButtonDownFcn', @StartWindowLevel);

	pData = getappdata(hFig, 'pData');
	set(hFig, 'Pointer', 'custom', 'PointerShapeCData', pData.windowLevel, 'PointerShapeHotSpot', [8 8])

	setappdata(hFig, 'strMode', 'Window/Level');
	UpdateStatusText(hFig);
end

% ------------------------------------------------------------------------------
function ToggleRoiDrawMode(hFig)
	if ~isempty(get(hFig, 'WindowButtonDownFcn')) && strcmp(func2str(get(hFig, 'WindowButtonDownFcn')), func2str(@StartRoiDraw))
		% In "Draw ROI" mode with no ROIs currently being drawn
		SetNormalMode(hFig);
		return
	elseif ~isempty(get(hFig, 'WindowButtonDownFcn')) && strcmp(func2str(get(hFig, 'WindowButtonDownFcn')), func2str(@InProgressRoiClick))
		% In "Draw ROI" mode, but an ROI is currently being drawn
		set(getappdata(hFig, 'hStatusText'), 'String', 'Error: Cannot exit ROI draw mode while drawing an ROI', 'Color', 'r');
		return
	end

	SetNormalMode(hFig);  % Clean up any cruft from other modes (e.g. display from circle tool)
	set(hFig, 'WindowButtonDownFcn', @StartRoiDraw);

	setptr(hFig, 'datacursor');

	setappdata(hFig, 'strMode', 'Draw New ROI');
	UpdateStatusText(hFig);
end

% ------------------------------------------------------------------------------
function ToggleColormap(hFig, bCustomMap)

	if verLessThan('matlab','8.4.0')
		currentColormap = get(hFig, 'Colormap');
	else
		CurrentAxes     = get(hFig, 'CurrentAxes');
		currentColormap = colormap(CurrentAxes);
	end
	prefs = getappdata(hFig, 'prefs');
	hAxes = reshape(getappdata(hFig, 'hAxes'), 1, []);
	hAxes(hAxes == 0) = [];
	
	% Default colormap if we can't identify current colormap or input fails
	newCMap = gray(size(currentColormap,1));

	if exist('bCustomMap', 'var') && ~isempty(bCustomMap) && bCustomMap
		% Custom colormap using user defined function
		strCMap = input('Function name for custom colormap: ', 's');
		if exist(strCMap, 'file')
			try
				newCMap = eval(sprintf('%s(%d)', strCMap, size(currentColormap,1)));
			catch ME
				% Lazy, so try to support this
				if strcmp(ME.identifier, 'MATLAB:dispatcher:InexactCaseMatch')
					[~, strCMap] = fileparts(which(strCMap));
					newCMap = eval(sprintf('%s(%d)', strCMap, size(currentColormap,1)));
				else
					rethrow(ME)
				end
			end
		else
			disp(['Could not find function: ' strCMap]);
		end
	else
		% Identify current colormap
		bMatchColormap = cellfun(@(x) all(all((x(size(currentColormap,1)) - currentColormap) == 0)), prefs.colormaps);
		indCurrColormap = find(bMatchColormap,1);

		if ~isempty(indCurrColormap)
			indNextColormap = mod(indCurrColormap, numel(prefs.colormaps))+1;
			newCMap = prefs.colormaps{indNextColormap}(size(currentColormap,1));
		end
	end

	if verLessThan('matlab','8.4.0')
		set(hFig, 'Colormap', newCMap);
		colormap(hFig, newCMap);
	else
		arrayfun(@(x) colormap(x,newCMap), hAxes);
	end

	% Apply to m-mode image
	roiPlot = getappdata(hFig, 'roiPlot');
	if ~isempty(roiPlot) && ~isempty(roiPlot.hFig)
		set(roiPlot.hFig, 'Colormap', newCMap);
	end
end

% ------------------------------------------------------------------------------
function ToggleStaticRoi(hFig)
	staticRoi = getappdata(hFig, 'staticRoi');
	if staticRoi
		setappdata(hFig, 'staticRoi', false);
	else
		setappdata(hFig, 'staticRoi', true);
	end

	UpdateStatusText(hFig)
end

% ------------------------------------------------------------------------------
function SetNormalMode(hFig)
	CurrentAxes = get(hFig, 'CurrentAxes');

	% If applicable, delete preview circle for Circle tool
	hCircle = getappdata(CurrentAxes, 'hCircle');
	if ~isempty(hCircle)
		delete(hCircle);
	end
	setappdata(CurrentAxes, 'hCircle', []);

	% Clear the undo ROIs for the nudge tool
	setappdata(CurrentAxes, 'prevRoi', [])
	
	set(hFig, 'WindowButtonDownFcn',   @MouseDownSelect);
	set(hFig, 'WindowButtonUpFcn',     '');
	set(hFig, 'WindowButtonMotionFcn', '');
	set(hFig, 'Pointer',               'arrow');
	setappdata(hFig, 'strMode', 'Normal');
	UpdateStatusText(hFig);
end

%% -----------------------------------------------------------------------------
function CreateInfoFigure(hFig)

	try
		hAxes = getappdata(hFig, 'hAxes');
		extras = getappdata(hFig, 'extras');
		CurrentAxes = get(hFig, 'CurrentAxes');
		CurrentSubplot = (hAxes == CurrentAxes);

		currentTime = getappdata(CurrentAxes, 'currentTime');
		currentSlice = getappdata(CurrentAxes, 'currentSlice');		

		info = extras.headers{CurrentSubplot}{currentSlice}{currentTime};

		strSCD = ShowCommonDicom(info);

		if isfield(info, 'ascconv')
			% Do nothing, it's already been extracted
		elseif isfield(info, 'Private_0029_1020')
			ascconv = ExtractSiemensHeader(char(info.Private_0029_1020)');
		end
	
		outString = [evalc('info') ...
		             sprintf('\n') evalc('ascconv') ...
		             sprintf('\n') strSCD];

		hFigInfo = figure;
		hText = uicontrol(hFigInfo, ...
											'style', 'list', ...
											'string', outString, ...
											'FontName', 'Lucida Console', ...
											'Units', 'normalized', 'Position', [0 0 1 1], ...
											'Min', 0, ...
											'Max', sum(outString == sprintf('\n'))+1);

		extent = get(hText, 'Extent');
		figPos = get(hFigInfo, 'Position');
		figPos(3) = figPos(3) * extent(3);
		set(hFigInfo, 'Position', figPos);

		% Scroll to end using a little java hackery
		jhEdit = findjobj(hText);
		jEdit = jhEdit.getComponent(0).getComponent(0);
		nLines = jEdit.getModel.getSize;
		jEdit.ensureIndexIsVisible( nLines - 1 );

		% And make it easily dismissible
		set(hFigInfo, 'KeyPressFcn',         @CloseInfoFig);
		set(hText,    'KeyPressFcn',         @CloseInfoFig);
	catch ME
		rethrow(ME)
	end
end

function CloseInfoFig(src, evt)
	if strcmp(evt.Key, 'escape')
		if strcmp(get(src, 'Type'), 'figure')
			close(src);
		else
			close(get(src, 'Parent'));
		end
	end
end

function ascconv = ExtractSiemensHeader(str)
	% The Siemens ASCII header is stored in the DICOM header under Private_0029_1020
	% and its beginning and end are marked by "### ASCCONV BEGIN ###" and 
	% "### ASCCONV END ###" respectively.  However:
	% - There are occasionally corrupt header fragments, that have an
	%   "### ASCCONV BEGIN ###" without a matching ### ASCCONV END ###"
	% - There may be # characters within the header to denote comments.
	%
	% A simple regular expression can't handle this.  The following approach uses
	% non-greedy matching to strip out header fragments that have a beginning
	% (with no end), but does not work for fragments with an end (and no beginning)
% 	ascConvText = regexp(str(end:-1:1), '### DNE VNOCCSA ###\n(?<vnoccsa>.+?)\n###\s+NIGEB VNOCCSA ###', 'once', 'names');
	if iscolumn(str)
		str = str';
	end
	ascConvText = regexp(str(end:-1:1), '### DNE VNOCCSA ###\n(?<vnoccsa>.+?)\n###[^\n]+NIGEB VNOCCSA ###', 'once', 'names');
	if ~isempty(ascConvText)
		ascconv = ascConvText.vnoccsa(end:-1:1);
	else
		ascconv = '';
	end
	
	% Add an extra check for an ASCCONV END without an ASCCONV BEGIN
	indEnd = findstr(ascconv, '### ASCCONV END ###');
	if ~isempty(indEnd)
		ascconv = ascconv(1:indEnd-2);
	end
end

%%
function TogglePlotIntersect(hFig)
	prefs       = getappdata(hFig, 'prefs');
	extras      = getappdata(hFig, 'extras');
	hAxes       = getappdata(hFig, 'hAxes');
	CurrentAxes = get(       hFig, 'CurrentAxes');

	% Check to see if there is DICOM information
	if isempty(extras) || ~isfield(extras, 'headers')
		return
	end

	% Abort if only one series displayed
	if (numel(hAxes) == 1)
		return
	end
	
	hIntPlot = getappdata(CurrentAxes, 'hIntPlot');
	if ~isempty(hIntPlot)
		% Clean up old plots
		for iSub = 1:numel(hAxes)
			if (hAxes(iSub) == 0)
				continue
			end
			hIntPlot = getappdata(hAxes(iSub), 'hIntPlot');
% 			delete(hIntPlot((1:numel(hAxes)) ~= iSub))
			delete(hIntPlot(hIntPlot ~= 0))
			setappdata(hAxes(iSub), 'hIntPlot', []);
			
			% Reset text color
			hText = getappdata(hAxes(iSub), 'hText');
			set(hText, 'Color', 'w')
% 			set(hText, 'EdgeColor', 'none')
		end
		return
	end

	% Collect the slice position information from all slices
	infoPos = repmat(struct('ImageOrientationPatient', [], 'ImagePositionPatient', [], 'Rows', [], 'Columns', [], 'PixelSpacing', []), [1 numel(hAxes)]);
	for iSub = 1:numel(hAxes)
		if (hAxes(iSub) == 0)
			continue
		end
		currentSlice = getappdata(hAxes(iSub), 'currentSlice');
		currentTime  = getappdata(hAxes(iSub), 'currentTime');

		try
			infoPos(iSub).ImageOrientationPatient = extras.headers{iSub}{currentSlice}{currentTime}.ImageOrientationPatient;
			infoPos(iSub).ImagePositionPatient    = extras.headers{iSub}{currentSlice}{currentTime}.ImagePositionPatient;
			infoPos(iSub).Rows                    = extras.headers{iSub}{currentSlice}{currentTime}.Rows;
			infoPos(iSub).Columns                 = extras.headers{iSub}{currentSlice}{currentTime}.Columns;
			infoPos(iSub).PixelSpacing            = extras.headers{iSub}{currentSlice}{currentTime}.PixelSpacing;
		catch
			% Maybe slice and time have been switched
			infoPos(iSub).ImageOrientationPatient = extras.headers{iSub}{currentTime}{currentSlice}.ImageOrientationPatient;
			infoPos(iSub).ImagePositionPatient    = extras.headers{iSub}{currentTime}{currentSlice}.ImagePositionPatient;
			infoPos(iSub).Rows                    = extras.headers{iSub}{currentTime}{currentSlice}.Rows;
			infoPos(iSub).Columns                 = extras.headers{iSub}{currentTime}{currentSlice}.Columns;
			infoPos(iSub).PixelSpacing            = extras.headers{iSub}{currentTime}{currentSlice}.PixelSpacing;
		end
	end

	% Plot intersections
	for iSub = 1:numel(hAxes)
		if (hAxes(iSub) == 0)
			continue
		end

		clear hIntPlot
		for jSub = 1:numel(hAxes)
			if (hAxes(jSub) == 0) || (iSub == jSub)
				continue
			end
			[xVert, yVert] = GetLpsIntersectFromDicom(infoPos(jSub), infoPos(iSub));

			hIntPlot(jSub) = plot(hAxes(iSub), xVert, yVert, 'Color', prefs.roiColors(mod(jSub-1, size(prefs.roiColors,1))+1,:));

		end
		setappdata(hAxes(iSub), 'hIntPlot', hIntPlot)
		
		% Change colour of text to identify which line is which subplot
		hText = getappdata(hAxes(iSub), 'hText');
		set(hText, 'Color', prefs.roiColors(mod(iSub-1, size(prefs.roiColors,1))+1,:))
% 		set(hText, 'EdgeColor', prefs.roiColors(iSub,:))
% 		set(hText, 'LineWidth', 2)
	end
end

%%
function UpdatePlotIntersect(hFig)
% 	prefs       = getappdata(hFig, 'prefs');
	extras      = getappdata(hFig, 'extras');
	hAxes       = getappdata(hFig, 'hAxes');
	CurrentAxes = get(       hFig, 'CurrentAxes');

	hIntPlot = getappdata(CurrentAxes, 'hIntPlot');
	if isempty(hIntPlot)
		return
	end

	% Collect the slice position information from all slices
	infoPos = repmat(struct('ImageOrientationPatient', [], 'ImagePositionPatient', [], 'Rows', [], 'Columns', [], 'PixelSpacing', []), [1 numel(hAxes)]);
	for iSub = 1:numel(hAxes)
		if (hAxes(iSub) == 0)
			continue
		end
		currentSlice = getappdata(hAxes(iSub), 'currentSlice');
		currentTime  = getappdata(hAxes(iSub), 'currentTime');

		try
			infoPos(iSub).ImageOrientationPatient = extras.headers{iSub}{currentSlice}{currentTime}.ImageOrientationPatient;
			infoPos(iSub).ImagePositionPatient    = extras.headers{iSub}{currentSlice}{currentTime}.ImagePositionPatient;
			infoPos(iSub).Rows                    = extras.headers{iSub}{currentSlice}{currentTime}.Rows;
			infoPos(iSub).Columns                 = extras.headers{iSub}{currentSlice}{currentTime}.Columns;
			infoPos(iSub).PixelSpacing            = extras.headers{iSub}{currentSlice}{currentTime}.PixelSpacing;
		catch
			% Maybe slice and time have been switched
			infoPos(iSub).ImageOrientationPatient = extras.headers{iSub}{currentTime}{currentSlice}.ImageOrientationPatient;
			infoPos(iSub).ImagePositionPatient    = extras.headers{iSub}{currentTime}{currentSlice}.ImagePositionPatient;
			infoPos(iSub).Rows                    = extras.headers{iSub}{currentTime}{currentSlice}.Rows;
			infoPos(iSub).Columns                 = extras.headers{iSub}{currentTime}{currentSlice}.Columns;
			infoPos(iSub).PixelSpacing            = extras.headers{iSub}{currentTime}{currentSlice}.PixelSpacing;
		end
	end

	% Update intersection plots
	for iSub = 1:numel(hAxes)
		if (hAxes(iSub) == 0)
			continue
		end

		hIntPlot = getappdata(hAxes(iSub), 'hIntPlot');
		for jSub = 1:numel(hAxes)
			if (hAxes(jSub) == 0) || (iSub == jSub)
				continue
			end
			[xVert, yVert] = GetLpsIntersectFromDicom(infoPos(jSub), infoPos(iSub));

			set(hIntPlot(jSub), 'XData', xVert, 'YData', yVert)
		end
	end
end


function [xVert, yVert, lVert1, lVert2] = GetLpsIntersectFromDicom(infoOther, infoRef)
% GetLpsIntersectFromDicom  Calculate interception between DICOM images
%
% Note: Adapted from PlotLpsFromDicom

	xVert = nan(length(infoOther), 2);
	yVert = nan(length(infoOther), 2);

	if ~iscolumn(  infoRef.ImageOrientationPatient),   infoRef.ImageOrientationPatient =   infoRef.ImageOrientationPatient'; end
	if ~iscolumn(infoOther.ImageOrientationPatient), infoOther.ImageOrientationPatient = infoOther.ImageOrientationPatient'; end
	if ~iscolumn(  infoRef.ImagePositionPatient   ),   infoRef.ImagePositionPatient    =   infoRef.ImagePositionPatient';    end
	if ~iscolumn(infoOther.ImagePositionPatient   ), infoOther.ImagePositionPatient    = infoOther.ImagePositionPatient';    end
	
	%  A transformation matrix transRef is constructed:
	%   [ i1 j1 k1 t1 ]     i is the x direction vector
	%   [ i2 j2 k2 t2 ]     j is the y direction vector
	%   [ i3 j3 k3 t3 ]     k is the z direction vector
	%   [ 0  0  0  1  ]     t is the origin
	normalRef = cross(infoRef.ImageOrientationPatient(1:3), infoRef.ImageOrientationPatient(4:6));
	transRef = [[infoRef.ImageOrientationPatient(1:3), infoRef.ImageOrientationPatient(4:6), normalRef, infoRef.ImagePositionPatient]; 0 0 0 1];

	for i = 1:length(infoOther)
		% Convert the normal and origin for the other images into the logical
		% coordinates for the reference image by multiplying by the inverse of 
		% the transformation matrix
		normalOther   = cross(infoOther(i).ImageOrientationPatient(1:3), infoOther(i).ImageOrientationPatient(4:6));
		normalOther_Ref = transRef(1:3,1:3) \ normalOther;
		originOther_Ref = transRef          \ [infoOther(i).ImagePositionPatient; 1];

		% The reference image can be represented by the equation
		%   ax + by + cz + d = 0
		a = normalOther_Ref(1);
		b = normalOther_Ref(2);
		c = normalOther_Ref(3);
		d = -(a*originOther_Ref(1) + b*originOther_Ref(2) + c*originOther_Ref(3));

		% Check for parallel images
		if (a == 0) && (b == 0)
			continue
		end

		% Set z = 0 and calculate the intersection line with either
		%   y = -1/b * (a*x + d)
		%   x = -1/a * (b*y + d)
		if (a == 0)
			xLine = double([1 infoRef.Columns]) * infoRef.PixelSpacing(1);
			yLine = -1/b * (a*xLine + d);
		else
			yLine = double([1 infoRef.Rows]) * infoRef.PixelSpacing(2);
			xLine = -1/a * (b*yLine + d);
		end

		% Save vertices
		xVert(i,:) = xLine / infoRef.PixelSpacing(1);
		yVert(i,:) = yLine / infoRef.PixelSpacing(2);

		lVert1 = transRef * [xVert(i,1) yVert(i,1) 0 1]';
		lVert2 = transRef * [xVert(i,2) yVert(i,2) 0 1]';
	end
end

%% -----------------------------------------------------------------------------
function CloseSimpleViewer(src,evt)
	prefs = getappdata(src, 'prefs');

	if prefs.confirmClose
		in = questdlg('Are you sure you want to close this window?', 'SimpleViewer', 'Yes', 'No', 'No');
		if ~strcmp(in, 'Yes')
			return
		end
	end

	roiPlot = getappdata(src, 'roiPlot');
	if ishandle(roiPlot.hFig)
		delete(roiPlot.hFig);
	end

	% Apply to SimpleViewer3D, if applicable
	hSV3D = getappdata(src, 'hSV3D');
	if ~isempty(hSV3D) && ishandle(hSV3D)
		close(hSV3D)
	end

	if ishandle(src)
		delete(src);
	end
end

% ------------------------------------------------------------------------------
function UpdateStatusText(hFig)
% Updates status text
	staticRoi   = getappdata(hFig, 'staticRoi');
	strMode     = getappdata(hFig, 'strMode');
	strError    = getappdata(hFig, 'strError');
	hStatusText = getappdata(hFig, 'hStatusText');
	prefs       = getappdata(hFig, 'prefs');

	% Build status text string
	strStatus = '';
	if ~isempty(strError)
		strStatus = strError;
	end

	if prefs.showStaticRoiStatus && ~staticRoi
		if isempty(strStatus)
			strStatus = 'StaticRoi is OFF';
		else
			strStatus = [strStatus sprintf('\nStaticRoi is OFF')];
		end
	end

	if ~strcmp(strMode, 'Normal')
		if isempty(strStatus)
			strStatus = ['Mode: ' strMode];
		else
			strStatus = [strStatus sprintf('\nMode: ') strMode];
		end
	end

	% Get status text colour
	if ~isempty(strError)
		statusColor = 'r';
	elseif ~strcmp(strMode, 'Normal')
		statusColor = 'y';
	else
		statusColor = 'w';
	end

	if isempty(strStatus)
		set(hStatusText, 'Visible', 'off')
	else
		set(hStatusText, 'Position', [0.01 0.01 0], ...
										 'Color',    statusColor, ...
										 'String',   strStatus, ...
										 'Visible',  'on');
	end
end

% ------------------------------------------------------------------------------
function lims = AutoLims(data, tol)
% Calculates limits to clip 2.5% of data at both ends

	if ~exist('tol', 'var')
% 		tol = [0.01 0.99];
		tol = [0.025 0.975];
	end

	sortedData = sort(data(~isnan(data) & ~isinf(data)));
	
	if isempty(sortedData)
		lims = nan(1,2);
		return
	end

	inds = [max([1                 round(numel(sortedData)*tol(1))]), ...
	        min([numel(sortedData) round(numel(sortedData)*tol(end))])];

	lims = sortedData(inds);

	% CLim must be increasing
	if numel(unique(lims)) == 1
		lims(end) = lims(end) + 1e-9;
	end
end

% ------------------------------------------------------------------------------
function scale = CalculateCurrentScale(hFig)
% Computes the scale for the current axis's image

	prefs = getappdata(hFig, 'prefs');

	hCurrentAxis = get(hFig, 'CurrentAxes');
	set(hCurrentAxis, 'Units', 'Pixels')

	position = get(hCurrentAxis, 'Position');
	imgSize  = size(get(findobj(hCurrentAxis, 'Type', 'image'), 'CData'));

	% Use this to zoom based on the displayed zoom, not the full image
	imgSize = round([diff(get(hCurrentAxis, 'YLim')) diff(get(hCurrentAxis, 'XLim'))]);

	set(hCurrentAxis, 'Units', 'normalized')
	
	scaleX = position(3) / imgSize(2) / (prefs.PixelSpacing(2) / prefs.PixelSpacing(1));
	scaleY = position(4) / imgSize(1);

% 	[scaleX scaleY]
	scale = min([scaleX scaleY]);
end

% ------------------------------------------------------------------------------
function SetFigureScale(hFig, scale, minSize)
% Sets a scale factor for the current image and applies it to the figure

	prefs = getappdata(hFig, 'prefs');
	CurrentAxes = get(hFig, 'CurrentAxes');
	szImg  = size(get(getappdata(CurrentAxes, 'hImage'), 'CData'));
	szAxes = size(getappdata(hFig, 'hAxes'));
	
	% Use this to zoom based on the displayed zoom, not the full image
	szImg = round([diff(get(CurrentAxes, 'YLim')) diff(get(CurrentAxes, 'XLim'))]);
	
	% Abort if the current frame has no image data
	if all(szImg == 0)
		return
	end

	% Remember, hAxes is the transpose of the subplot layout
	nCols = szAxes(1);
	nRows = szAxes(2);

	newHeight = szImg(1) * scale * nRows;
	newWidth  = szImg(2) * scale * nCols;

	% Get screen size
	previousUnits = get(0, 'Units');
	set(0, 'Units', 'pixels');
	screenPos = get(0, 'ScreenSize');  % [left bottom width height]

	% Mac has a non-hideable file menu
	if ismac
% 		screenPos(4) = screenPos(4) - 45;
		screenPos(4) = screenPos(4) - 71;
	end

	% Enforce a minimum size
	if exist('minSize', 'var') && ~isempty(minSize)
		ppi = get(0, 'ScreenPixelsPerInch');
		if (szImg(1) * nRows * scale)/ppi < minSize
			scale = (minSize*ppi) / (szImg(1)*nRows);
		end

		if (szImg(2) * nCols * scale * prefs.PixelSpacing(2) / prefs.PixelSpacing(1))/ppi < minSize
			scale = (minSize*ppi*prefs.PixelSpacing(1)/prefs.PixelSpacing(2)) / (szImg(2)*nCols);
		end
	end
	newHeight = szImg(1) * scale * nRows;
	newWidth  = szImg(2) * scale * nCols;
	
	% Reduce scale to fit screen dimensions
	scaleY = scale;
	if newHeight > screenPos(4)
		scaleY = screenPos(4) / (szImg(1) * nRows);
	end

	scaleX = scale;
	if newWidth > screenPos(3)
		scaleX = screenPos(3) / (szImg(2) * nCols * prefs.PixelSpacing(2) / prefs.PixelSpacing(1));
	end

% 	[scaleX scaleY]
	scale     = min([scaleX scaleY]);
	newHeight = szImg(1) * scale * nRows;
	newWidth  = szImg(2) * scale * nCols * prefs.PixelSpacing(2) / prefs.PixelSpacing(1);

	% Try to keep the top left fixed, but if it would push the figure out of screen, move it
	windowPos = get(hFig, 'Position');

	if (windowPos(1) + newWidth) > screenPos(3)
		windowPos(1) = screenPos(3) - newWidth;
	end

	bottom = windowPos(2) + windowPos(4) - newHeight;
	if (bottom < 0)
		windowPos(2) = 0;
	elseif (bottom + newHeight) > screenPos(4)
		windowPos(2) = screenPos(4) - newHeight;
	else
		windowPos(2) = bottom;
	end

	windowPos(3) = newWidth;
	windowPos(4) = newHeight;

	% Set window position
	set(hFig, 'Position', windowPos)

	% Reset units
	set(0, 'Units', previousUnits);
end

% ------------------------------------------------------------------------------
function strText = ComposeText(strTitle, iSlice, nSlices, iTime, nTimes)
% Compose a text string to display for each image, omitting portions that are not useful

	strSlice = '';
	if nSlices > 1
		strSlice = sprintf('Slice %d/%d', iSlice, nSlices);
	end

	strTime = '';
	if nTimes > 1
		strTime = sprintf('Time %d/%d', iTime, nTimes);
	end

	% The second case covers the "both empty" case as well
	if ~isempty(strSlice) && ~isempty(strTime)
		strSliceTime = [strSlice ', ' strTime];
	elseif isempty(strSlice)
		strSliceTime = strTime(6:end);
	elseif isempty(strTime)
		strSliceTime = strSlice(7:end);
	end

	% Add the title if it exists
	if ~isempty(strTitle) && ~isempty(strSliceTime)
		strText = [strTitle sprintf('\n') strSliceTime];
	elseif isempty(strTitle)
		strText = strSliceTime;
	elseif isempty(strSliceTime)
		strText = strTitle;
	end
end

% ------------------------------------------------------------------------------
function aspectRatio = GetAspectRatio()
% Get the current screen's aspect ratio
	previousUnits = get(0, 'Units');
	set(0, 'Units', 'pixels');
	screenPos = get(0, 'ScreenSize');  % [left bottom width height]
	set(0, 'Units', previousUnits);

	aspectRatio = screenPos(3) / screenPos(4);

	% It can occasionally be screwed up with multi-monitor setups
	if aspectRatio > (16/9)
		aspectRatio = 16/9;
	end
end

% ------------------------------------------------------------------------------
function CreateCustomPointers(hFig)
	pData.panZoom = ...
    [NaN NaN NaN NaN   1   1   1   1 NaN NaN NaN NaN   1 NaN NaN NaN; ...
     NaN NaN   1   1 NaN   2 NaN   2   1   1 NaN   1   2   1 NaN NaN; ...
     NaN   1   2 NaN   2 NaN   2 NaN   2 NaN   1   1   2   1   1 NaN; ...
     NaN   1 NaN   2 NaN   2 NaN   2 NaN   1   2   2   2   2   2   1; ...
       1 NaN   2 NaN   2 NaN   2 NaN   2 NaN   1   1   2   1   1 NaN; ...
       1   2 NaN   2 NaN   2 NaN   2 NaN   2 NaN   1   2   1 NaN NaN; ...
       1 NaN   2 NaN   2 NaN   2 NaN   2 NaN   2   1   1 NaN NaN NaN; ...
       1   2 NaN   2 NaN   2 NaN   2 NaN   2 NaN   1 NaN NaN NaN NaN; ...
     NaN   1   2 NaN   2 NaN   2 NaN   2 NaN   1 NaN NaN NaN NaN NaN; ...
     NaN   1 NaN   2 NaN   2 NaN   2 NaN   2   1   2 NaN NaN NaN NaN; ...
     NaN NaN   1   1   2 NaN   2 NaN   1   1   1   1   2 NaN NaN NaN; ...
     NaN NaN NaN NaN   1   1   1   1 NaN   2   1   1   1   2 NaN NaN; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   2 NaN; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   2; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   2];

	pData.panZoom = ...
    [NaN NaN NaN NaN   1   1   1   1 NaN NaN NaN NaN   1 NaN NaN NaN; ...
     NaN NaN   1   1 NaN   2 NaN   2   1   1 NaN   1   2   1 NaN NaN; ...
     NaN   1   2 NaN   2 NaN   2 NaN   2 NaN   1   1   2   1   1 NaN; ...
     NaN   1 NaN   2 NaN   2 NaN   2 NaN   1   2   2   2   2   2   1; ...
       1 NaN   2 NaN   2 NaN   2 NaN   2 NaN   1   1   2   1   1 NaN; ...
       1   2 NaN   2 NaN   2 NaN   2 NaN   2 NaN   1   2   1 NaN NaN; ...
       1 NaN   2 NaN   2 NaN   2 NaN   2 NaN   2   1   1 NaN NaN NaN; ...
       1   2 NaN   2 NaN   2 NaN   2 NaN   2 NaN   1 NaN NaN NaN NaN; ...
     NaN   1   2 NaN   2 NaN   2 NaN   2 NaN   1 NaN NaN NaN NaN NaN; ...
     NaN   1 NaN   2 NaN   2 NaN   2 NaN   2   1   2 NaN NaN NaN NaN; ...
     NaN NaN   1   1   2 NaN   2 NaN   1   1   1   1   2 NaN NaN NaN; ...
     NaN NaN NaN NaN   1   1   1   1 NaN   2   1   1   1   2 NaN NaN; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   2 NaN; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   2; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1; ...
     NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   2   1   2];
	% Quirky... the 'glass' pointer is shifted down and right by 1 pixel from what it claims.
	% Do this so that they line up.
	pData.panZoom = [nan(1,16); pData.panZoom(1:end-1,:)];
	pData.panZoom = [nan(16,1) pData.panZoom(:,1:end-1)];
	
	pData.windowLevel = ...
    [NaN NaN NaN NaN NaN   1   1   1   2   2   2 NaN NaN NaN NaN NaN; ...
     NaN NaN NaN   1   1   2   2   2   1   1   1   2   2 NaN NaN NaN; ...
     NaN NaN   1   2   2   2   2   2   1   1   1   1   1   2 NaN NaN; ...
     NaN   1   2   2   2   2   2   2   1   1   1   1   1   1   2 NaN; ...
     NaN   1   2   2   2   2   2   2   1   1   1   1   1   1   2 NaN; ...
       1   2   2   2   2   2   2   2   1   1   1   1   1   1   1   2; ...
       1   2   2   2   2   2   2   2   1   1   1   1   1   1   1   2; ...
       1   2   2   2   2   2   2   2   1   1   1   1   1   1   1   2; ...
       1   2   2   2   2   2   2   2   1   1   1   1   1   1   1   2; ...
       1   2   2   2   2   2   2   2   1   1   1   1   1   1   1   2; ...
       1   2   2   2   2   2   2   2   1   1   1   1   1   1   1   2; ...
     NaN   1   2   2   2   2   2   2   1   1   1   1   1   1   2 NaN; ...
     NaN   1   2   2   2   2   2   2   1   1   1   1   1   1   2 NaN; ...
     NaN NaN   1   2   2   2   2   2   1   1   1   1   1   2 NaN NaN; ...
     NaN NaN NaN   1   1   2   2   2   1   1   1   2   2 NaN NaN NaN; ...
     NaN NaN NaN NaN NaN   1   1   1   2   2   2 NaN NaN NaN NaN NaN];
 
	pData.forceField = nan(16,16);
	pData.forceField([71 72 73 74 86 91 101 108 117 124 133 140 149 156 166 171 183 184 185 186]) = 1;
	pData.forceField([6 7 8 9 10 11 20 21 28 29 35 45 46 50 63 66 79 81 96 97 112 113 128 129 144 145 160 161 176 178 191 194 207 211 221 222 228 229 236 237 246 247 248 249 250 251]) = 2;

	setappdata(hFig, 'pData', pData);
end

function interpRoi = CubicInterp(roi, nPts, method)
% Perform a cubic interpolation on ROIs

	tol = 1e-6; % distances smaller than this are considered negligble

	if size(roi,1) <= 1
		% No interpolation possible on a point or a line
		interpRoi = roi;
		return
	end
	
	% Return if it's two identical points (occurs when first drawing an ROI)
	if (size(roi,1) == 2) && sum(abs(diff(roi))) < tol
		interpRoi = roi;
		return
	end

	if ~exist('nPts', 'var') || isempty(nPts)
		nPts = 50;%min(50, size(roi,1)*5);
	end
	
	if ~exist('method', 'var') || isempty(method)
		method = 'spline';
	end

	% Don't close polygon for line segments
	if all(isnan(roi(end,:)))
		roi = roi(1:end-1,:);
		bIsLine = true;
	else
		bIsLine = false;
		% Close the polygon if necessary
		if (sum(abs(diff(roi([1 end],:)))) > tol)
			roi = roi([1:end 1],:);
		end
	end

	% Remove duplicated points
	isDuplicate = all(abs(diff(roi,1)) < tol, 2);
	if any(isDuplicate)
		roi(isDuplicate,:) = [];
	end
	
	% Interpolate x and y independently, based on cumulative length along the polygon's circumference
	diffDist = sqrt(diff(roi(:,1)).^2 + diff(roi(:,2)).^2);
	cumDist = cat(1, 0, cumsum(diffDist));

	% Spline interpolation
	interpRoi = zeros(nPts,2);
	interpRoi(:,1) = interp1(cumDist, roi(:,1), linspace(0,cumDist(end),nPts), method);
	interpRoi(:,2) = interp1(cumDist, roi(:,2), linspace(0,cumDist(end),nPts), method);

	% This enforces a zero slope at the end (not necessarily better)
% 	interpRoi(:,1) = spline(cumDist, [0; roi(:,1); 0], linspace(0,cumDist(end),nPts));
% 	interpRoi(:,2) = spline(cumDist, [0; roi(:,2); 0], linspace(0,cumDist(end),nPts));

% 	% Try using polar coordinates instead
% 	com = mean(roi);
% 	
% 	rad = sqrt(sum(bsxfun(@minus, roi, com).^2,2));
% 	tht = unwrap(atan2((roi(:,2)-com(2)), roi(:,1)-com(1)));
% 
% 	radi = interp1(cumDist, rad, linspace(0,cumDist(end),nPts), method);
% 	thti = interp1(cumDist, tht, linspace(0,cumDist(end),nPts), method);
% 	interpRoi = bsxfun(@plus, [radi.*cos(thti); radi.*sin(thti)]', com);
% 	
% 	% Some other crazy way of interpolating (from the spline docs)
% 	interpRoi = spline(tht, roi', linspace(-pi, pi, nPts))';

	% Re-add NaN marker for line segments
	if bIsLine
		interpRoi = cat(1, interpRoi, [nan nan]);
	end
end

function HideExtras(hFig)
% Prepare figure for export by hiding text and changing the renderer
%
% 	if strcmp(get(getappdata(hFig, 'hAxesBox'), 'Visible'), 'on')
	if strcmp(get(hFig, 'Renderer'), 'zbuffer') || strcmp(get(hFig, 'Renderer'), 'opengl') % This is a more reliable indicator if everything has been hidden or not
		set(getappdata(hFig, 'hAxesBox'),    'Visible', 'off')
		set(getappdata(hFig, 'hStatusText'), 'Visible', 'off')
		for hAxis = reshape(getappdata(hFig, 'hAxes'), 1, [])
			set(getappdata(hAxis, 'hText'), 'Visible', 'off')
		end
		set(hFig, 'Renderer', 'painters')
	else
		% Return to normal
		set(getappdata(hFig, 'hAxesBox'),    'Visible', 'on')
		set(getappdata(hFig, 'hStatusText'), 'Visible', 'on')
		for hAxis = reshape(getappdata(hFig, 'hAxes'), 1, [])
			set(getappdata(hAxis, 'hText'), 'Visible', 'on')
		end
		if verLessThan('matlab','8.4.0')
			set(hFig, 'Renderer', 'zbuffer')
		else
			set(hFig, 'Renderer', 'opengl')
		end
	end

	% Toggle the SD bars on the plot too
	prefs = getappdata(hFig, 'prefs');
	prefs.showStdPlot = ~prefs.showStdPlot;
	setappdata(hFig, 'prefs', prefs);
	
	roiPlot = getappdata(hFig, 'roiPlot');
	if ~isempty(roiPlot.hFig) && ishandle(roiPlot.hFig)
		UpdateRoiFigure(hFig);
	end
end

function SaveSvWorkspace(hFig)
% Save all data from a SimpleViewer window
	hAxes = reshape(getappdata(hFig, 'hAxes'), 1, []);

	% Remove empty axes
	hAxes(hAxes == 0) = [];
	
	img  = arrayfun(@(x) getappdata(x, 'imgData'), hAxes, 'UniformOutput', false);
	roi  = arrayfun(@(x) getappdata(x, 'roi'),     hAxes, 'UniformOutput', false);
	txt  = arrayfun(@(x) getappdata(x, 'txtData'), hAxes, 'UniformOutput', false);
	cmap = arrayfun(@(x) colormap(x),              hAxes, 'UniformOutput', false);

	% The ability to toggle ROI mask means that imgData could be masked
	nonMaskedImgData = arrayfun(@(x) getappdata(x, 'nonMaskedImgData'), hAxes, 'UniformOutput', false);
	if any(cellfun(@(x) ~isempty(x), nonMaskedImgData))
		img = nonMaskedImgData;
	end

	prefs = getappdata(hFig, 'prefs');

	if strcmp(get(hAxes(1), 'CLimMode'), 'manual')
% 		CLim = get(hAxes(1), 'CLim');
		CLim = arrayfun(@(x) get(x, 'CLim'), hAxes, 'UniformOutput', false);
	else
		CLim = [];
	end

	szAxes = size(getappdata(hFig, 'hAxes'));
	tileSize = szAxes([2 1]);

% 	cmap = colormap(hFig);
	
	XLim = arrayfun(@(x) get(x, 'XLim'), hAxes, 'UniformOutput', false);
	YLim = arrayfun(@(x) get(x, 'YLim'), hAxes, 'UniformOutput', false);

	% Save extra information associated with the figure
	extras = getappdata(hFig, 'extras');

	% Dirty workaround for the fact that uigetfile() doesn't allow you to specify start location
	persistent lastPath
	if exist('lastPath', 'var') && ~isempty(lastPath)
		currentDir = pwd;
		cd(lastPath)
		cu = onCleanup(@()cd(currentDir)); % Return to pwd when this function ends (even if crash/dbstop)
	end

	% Prompt to save file
	[fileName, pathName] = uiputfile({'*.svmat','SimpleViewer workspace'; '*.mat','MAT-files'}, 'Save SimpleViewer workspace');
	if ischar(fileName)
		save('-v7.3', fullfile(pathName, fileName), 'img', 'CLim', 'roi', 'txt', 'tileSize', 'prefs', 'cmap', 'XLim', 'YLim', 'extras')
		lastPath = pathName;
	end
end

function colorsRGB = ConvertColorsToRGB(colors)
	if ischar(colors)
		colorsRGB = zeros(numel(colors),3);
		for i = 1:numel(colors)
			switch colors(i)
				case 'y'
					colorsRGB(i,:) = [1 1 0];
				case 'm'
					colorsRGB(i,:) = [1 0 1];
				case 'c'
					colorsRGB(i,:) = [0 1 1];
				case 'r'
					colorsRGB(i,:) = [1 0 0];
				case 'g'
					colorsRGB(i,:) = [0 1 0];
				case 'b'
					colorsRGB(i,:) = [0 0 1];
				case 'w'
					colorsRGB(i,:) = [1 1 1];
				case 'k'
					colorsRGB(i,:) = [0 0 0];
			end
		end
	else
		colorsRGB = colors;
	end
end

%% Helper functions
% These functions are copies of helper functions, but are duplicated here so
% that SimpleViewer can have less dependencies.  These should be refreshed
% regularly so that they can incorporate updates to these functions.

function prefs = GeneratePrefs(defaultPrefs, inputPrefs)
% Generate prefs from a default set and merging in explicitly defined fields
%
% Output 'prefs' is the same as 'defaultPrefs', but each field in 'inputPrefs'
% that also exists in 'defaultPrefs' uses the value from 'inputPrefs' instead.
% This is commonly used when generating preferences, using default values when
% not explicitly defined, but overriding otherwise. Fields in 'inputPrefs' but
% not in 'defaultPrefs' will generate an error.
%
% Revisions:
%  1.0.1:  - Initial release
%
% Kelvin Chow (kelvinc@ualberta.ca)
% Revision: 1.0.1  Date: 22 November 2010

	prefs = defaultPrefs;

	if isempty(inputPrefs)
		return
	end

	fields = fieldnames(inputPrefs);
	for iField = 1:length(fields)
		if isfield(prefs, fields{iField})
%         disp(sprintf('  Found field: %s', fields{iField}));
%         disp(sprintf('Current Value: %s', prefs.(fields{iField})));
%         disp(sprintf('    New Value: %s', inputPrefs.(fields{iField})));
			prefs.(fields{iField}) = inputPrefs.(fields{iField});
		else
			disp(sprintf('[GeneratePrefs] Error: Invalid field: %s', fields{iField}));
		end
	end
end

% ------------------------------------------------------------------------------
function theAxis = subplottight(nrows, ncols, thisPlot, fillRatio, varargin)
% subplottight  Create axes in tightly tiled positions.
%
% Syntax:
%  [theAxis] = subplottight(nrows, ncols, thisPlot, fillRatio, ...)
%
% Inputs:
%  nrows     - number of rows
%  ncols     - number of columns
%  thisPlot  - index of desired subplot (row-major numbering; same as subplot)
%  fillRatio - fraction of the full subplot to fill (e.g. 0.85.  default: 1)
%  ...       - property-value pairs to be passed to subplot
%
%  Outputs:
%   thisAxis - handle of subplot created (optional)
%
% Notes:
%  - Tight spacing of subplots may hide subplot titles.
%
% Revisions:
%  1.0.3:  - Allow additional property-value pairs to be passed to subplot
%  1.0.2:  - Set figure units to be normalized to allow for better resizing
%  1.0.1:  - Initial release
%
% Kelvin Chow (kelvin.chow@ualberta.ca)
% Department of Biomedical Engineering
% University of Alberta
% Revision: 1.0.3  Date: 24 November 2014
%
% Copyright (c) 2014 Kelvin Chow
%
% Please do not distribute without the author's permission.

	if ~exist('fillRatio', 'var') || isempty(fillRatio)
		fillRatio = 1;
	end
	
  width  = 1 / ncols;
  height = 1 / nrows;

	gapX = width  * (1-fillRatio);
	gapY = height * (1-fillRatio);
	
  [iCol, iRow] = ind2sub([ncols nrows], thisPlot);
  
  left   = (iCol - 1)     * width;
  bottom = (nrows - iRow) * height;

  h = subplot('Position', [left + gapX, ...
	                         bottom + gapY, ...
	                         width - 2*gapX, ...
	                         height - 2*gapY], ...
	                         varargin{:});
	set(h, 'Units', 'normalized');
% 	set(h,'LooseInset',get(h,'TightInset'))
  if (nargout == 0)
    return
  end
  theAxis = h;
end

% ------------------------------------------------------------------------------
function b = WindowLevel(a, window, level)
% WindowLevel  WindowLevels a matrix to the range [0,1]
% 
% Syntax:
%  b = WindowLevel(a, [window, level])
%
% Description:
%  Inputs:
%   a      - matrix to be windowleveled (can be n-dimensional)
%   window - pre-set window (optional)
%   level  - pre-set level (optional)
%  Outputs:
%   b      - matrix a after windowleveling (range [0,1])
%
% Revisions:
%  1.0.1:  - Initial release
%
% Kelvin Chow (kelvinc@ualberta.ca)
% Revision: 1.0.1  Date: 27 November 2008

	if nargin < 2
		window = std(double(a(:)))*4;
		level  = mean(double(a(:)));
	end

	b = double(a) - level;
	b = b ./ window + 0.5;

	b(b<0) = 0;
	b(b>1) = 1;
end

% ------------------------------------------------------------------------------
function h = overobj2(varargin)
%OVEROBJ2 Get handle of object that the pointer is over.
%   H = OVEROBJ2 searches all objects in the PointerWindow
%   looking for one that is under the pointer. Returns first
%   object handle it finds under the pointer, or empty matrix.
%
%   H = OVEROBJ2(FINDOBJ_PROPS) searches all objects which are
%   descendants of the figure beneath the pointer and that are
%   returned by FINDOBJ with the specified arguments.
%
%   Example:
%       h = overobj2('type','axes');
%       h = overobj2('flat','visible','on');
%
%   See also OVEROBJ, FINDOBJ
%
% From: http://undocumentedmatlab.com/blog/undocumented-mouse-pointer-functions/

	% Ensure root units are pixels
	oldUnits = get(0,'units');
	set(0,'units','pixels');

	% Get the figure beneath the mouse pointer & mouse pointer pos
	if verLessThan('matlab','8.4.0')
		fig = get(0,'PointerWindow'); 
		p = get(0,'PointerLocation');
	else
		fig = get(groot, 'CurrentFigure_I');
		p = get(groot, 'PointerLocation_I');
	end
	set(0,'units',oldUnits);

	% Look for quick exit (if mouse pointer is not over any figure)
	if fig==0,  h=[]; return;  end

	% Compute figure offset of mouse pointer in pixels
	figPos = getpixelposition(fig);
	x = (p(1)-figPos(1));
	y = (p(2)-figPos(2));

	% Loop over all figure descendants
	c = findobj(get(fig,'Children'),varargin{:});
	for h = c',
		 % If descendant contains the mouse pointer position, exit
		 r = getpixelposition(h);
		 if (x>r(1)) && (x<r(1)+r(3)) && (y>r(2)) && (y<r(2)+r(4))
				return
		 end
	end
	h = [];
end

%% -----------------------------------------------------------------------------
function value = reinput(strQuery, defaultValue, stringMode, defaultFormat, allowedRange)
% reinput  Prompt with user input, with default values.  Default value is set
%          when input is empty
% 
% Syntax:
%  value = reinput(strQuery, defaultValue [, stringMode])
%
% Inputs:
%  strQuery      - text string for query
%  defaultValue  - default value for query
%  stringMode    - 's' for character string input (optional)
%                  'bool' for boolean input
%  defaultFormat - sprintf format to display defaultValue with
%  allowedRange  - allowable range of input values.  format: [min max]
%
% Outputs:
%  value         - user input (defaultValue if input was blank)
%
% Revisions:
%  1.0.7:  - Support 'bool' type
%  1.0.6:  - Escape special characters when defaultValue is a string
%  1.0.5:  - Undo changes made in 1.0.4
%  1.0.4:  - Default back to 'input.m' behaviour if defaultValue is missing or empty
%  1.0.3:  - Support an allowable range of values
%  1.0.2:  - Support a default sprintf format for default input values
%  1.0.1:  - Initial release
%
% Kelvin Chow (kelvin.chow@ualberta.ca)
% Department of Biomedical Engineering
% University of Alberta
% Revision: 1.0.7  Date: 4 August 2012
%
% Copyright (c) 2014 Kelvin Chow
%
% Please do not distribute without the author's permission.

	if ~exist('stringMode', 'var')
		stringMode = [];
	end

	if ~exist('defaultFormat', 'var')
		defaultFormat = '%d';
	end
	
	if ~exist('allowedRange', 'var')
		allowedRange = [nan nan];
	end
	
	if strcmp(stringMode, 'bool')
		defaultFormat = '%s';
		if exist('defaultValue', 'var')
			if (defaultValue == true)
				defaultValue = 'true';
			elseif (defaultValue == false)
				defaultValue = 'false';
			end
		end
	end

	while(1)
		if strcmp(stringMode, 's') || strcmp(stringMode, 'bool')
			if ~exist('defaultValue', 'var') || isempty(defaultValue)
				value = input(sprintf('%s ['''']: ', strQuery), 's');
% 				value = 'Gadovist';
			else
				% Must escape special characters before displaying
				escDefaultValue = defaultValue;
				escDefaultValue = regexprep(escDefaultValue, '\\', '\\\\');
				escDefaultValue = regexprep(escDefaultValue, '%', '%%');
				value = input(sprintf('%s [''%s'']: ', strQuery, escDefaultValue), 's');
			end
		else
			if ~exist('defaultValue', 'var') || isempty(defaultValue)
				value = input(sprintf('%s []: ', strQuery));
			else
				strDefault = sprintf([defaultFormat ' '], defaultValue);
				value = input(sprintf('%s [%s]: ', strQuery, strDefault(1:end-1)));
			end
		end

		if isempty(value) && exist('defaultValue', 'var') && ~isempty(defaultValue)
			value = defaultValue;
		end

		if strcmp(stringMode, 's')
			break
		elseif strcmp(stringMode, 'bool')
			if any(strcmp(value, {'true', 't'}))
				value = 1;
				break
			elseif any(strcmp(value, {'false', 'f'}))
				value = 0;
				break
			else
				disp('[reinput] Error: inputs must be ''true'', ''false'', ''t'', or ''f''');
			end
		else
			if any(value < allowedRange(1)) || any(value > allowedRange(2))
				disp(sprintf(sprintf('[reinput] Error: inputs must be within the range %s to %s', defaultFormat, defaultFormat), allowedRange(1), allowedRange(2)));
				continue
			else
				break
			end
		end
	end
end

%%
function img = dicomread_wrapper(varargin)
	try
		img = dicomread(varargin{:});
	catch err
		if strcmp(err.identifier,'MATLAB:subsassigndimmismatch')
			% This is likely due to overlays with different sizes -- use this code
			% which skips parsing of overlays
			img = dicomread_no_overlay(varargin{:});
		else
			rethrow(err);
		end
	end
end

%% -----------------------------------------------------------------------------
function indFirst = FindFirstNonEmptyCell(cData)
	% Find index of first non-empty cell in a cell array
	indFirst = find(cellfun(@(x) ~isempty(x), cData),1);
end
%% Use this to programmatically move the mouse
% import java.awt.Robot;
% mouse = Robot;
% mouse.mouseMove(0, 0);
% screenSize = get(0, 'screensize');
% for i = 1: screenSize(4)
%       mouse.mouseMove(i, i);
%       pause(0.00001);
% end

%%
% ScreenSize = get(0, 'ScreenSize');
% figPos = get(gcf, 'Position');
% 
% if figPos(1) > ScreenSize(3)
%     figPos(1) = ScreenSize(3)-figPos(3);
%     set(gcf, 'Position', figPos)
% end