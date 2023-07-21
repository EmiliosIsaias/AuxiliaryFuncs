function [eqFlag] = compareSoftVers(struct1, struct2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fldNms1 = fieldnames(struct1); fldNms2 = fieldnames(struct2);
fldNms1(contains(fldNms1, 'AnalysisTime')) = [];
fldNms2(contains(fldNms2, 'AnalysisTime')) = [];
if all(string(fldNms1) == string(fldNms2))
    eqFlag = all(arrayfun(@(fn) struct1.(fn) == struct2.(fn), string(fldNms1)));
end