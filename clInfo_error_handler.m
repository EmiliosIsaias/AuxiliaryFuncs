function [out] = clInfo_error_handler(S, varargin)
warning( S.identifier, S.message)
out = table();
end