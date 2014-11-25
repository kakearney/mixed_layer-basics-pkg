function fieldTable = ecopathinputinfo;
%ECOPATHINPUTINFO Table of dimesions and units for ecopath input fields
%
% fieldTable = ecopathinputinfo;
%
% This function returns a table with information on the dimensions and
% units used for each field in a Ewe input structure.
%
% Output variables:
%
%   fieldTable: n x 3 cell array of strings
%               column 1:   names of each input variable
%               column 2:   dimensions of variable (ngroup = number of
%                           functional groups in the model, ndet = number
%                           of non-live detritus groups, ngear = number of
%                           fishing fleets in the model)
%               column 3:   units of each variable (M = mass, A = area or
%                           volume, T = time)

% Copyright 2008 Kelly Kearney

fieldTable = {...
    'ngroup'        '1 x 1'             'no unit'   
    'nlive'         '1 x 1'             'no unit'
    'ngear'         '1 x 1'             'no unit'
    'name'          'ngroup x 1'        'string'
    'pp'            'ngroup x 1'        'no unit'
    'areafrac'      'ngroup x 1'        'no unit'
    'bh'            'ngroup x 1'        'M/A'
    'pb'            'ngroup x 1'        '1/T'
    'qb'            'ngroup x 1'        '1/T'
    'ee'            'ngroup x 1'        'no unit'
    'ge'            'ngroup x 1'        'no unit'
    'gs'            'ngroup x 1'        'no unit'
    'dtImp'         'ngroup x 1'        'M/A/T'
    'b'             'ngroup x 1'        'M/A'
    'dc'            'ngroup x ngroup'   'no unit'
    'df'            'ngroup x ndet'     'no unit'
    'immig'         'ngroup x 1'        'M/A/T'
    'emig'          'ngroup x 1'        'M/A/T'
    'emigRate'      'ngroup x 1'        '1/T'
    'ba'            'ngroup x 1'        'M/A/T'
    'baRate'        'ngroup x 1'        '1/T'
    'fleetname'     'ngear x 1'         'string'   
    'landing'       'ngroup x ngear'    'M/A/T'  
    'discard'       'ngroup x ngear'    'M/A/T' 
    'discardFate'   'ngear x ndet'      'no unit'
    'maxrelpb'      'ngroup x 1'        'no unit'
    'maxrelfeed'    'ngroup x 1'        'no unit'
    'feedadj'       'ngroup x 1'        'no unit' 
    'fracsens'      'ngroup x 1'        'no unit'    
    'predeffect'    'ngroup x 1'        'no unit'
    'densecatch'    'ngroup x 1'        'no unit'
    'qbmaxqb0'      'ngroup x 1'        'no unit'
    'switchpower'   'ngroup x 1'        'no unit'
    'kv'            'ngroup x ngroup'   'no unit'
    'x'             'ngroup x ngroup'   'no unit'
    'd'             'ngroup x ngroup'   'no unit'
    'theta'         'ngroup x ngroup'   'no unit'};

