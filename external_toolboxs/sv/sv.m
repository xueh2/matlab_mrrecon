function sv(varargin)
	% sv  Short alias to SimpleViewer
	if nargin > 0
		figure
	end
	SimpleViewer(varargin{:});
end