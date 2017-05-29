#pragma once
#include "StructDefinitions.h"

/**
 * Implementation of HSMwfwconvert.m.
 *@param model the HSM model acquired from training
 *@param audioFeature the feature of the input audio to be converted in place
 *@remark In the whole MATLAB program, the last parameter modo is always 'fwa'. Therefore, it is neglected here.
*/
void HSMwfwconvert(const HSMModel& model, PicosStructArray& audioFeature);



