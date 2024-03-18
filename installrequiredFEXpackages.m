% installrequiredFEXpackages -- a SAMPLE SCRIPT 
% illustrating the use of the 'requireFEXpackage' function
% 
% Install required FEX packages
% 
% Add this script and the function 'requireFEXpackage.m' to your toolbox, 
% modify the 'WhatIsRequired',
% and instruct the user of your toolbox to run this script once 
% before using your toolbox.
% 
% (C) Igor Podlubny 2011.


% List the functions written by other FEX authors 
% and the ID numbers of the corresponding FEX packages (submissions).
% Each row must contain the name of the function (string)
% and the ID of the FEX packages to which it belongs (integer number).
% The two rows below in 'WhatIsRequired' array 
% illustrate how this must be done. Edit it as you need.
WhatIsRequired = {
'fminsearchbnd', 8277; % function 'fminsearchbnd' from FEX submission 8277 is needed 
'mlf', 8738;           % function 'mlf' from FEX submission 8738 is needed, too
'routh', 58            % function 'routh' -- single m-file without BSD licence
};


% The subsequent code processes the list given in 'WhatIsRequired'

% Number of required packages
n= size(WhatIsRequired,1);

for k=1:n
   % read function name and FEX package ID
   Fun = WhatIsRequired(k,1);                       
   Package = WhatIsRequired(k,2);
   % if the function does not exist in your Matlab path,
   % downloaded and install the FEX package containing that function.
   if ~(exist(char(Fun), 'file') == 2)              
        P = requireFEXpackage(cell2mat(Package));   
   end
end






