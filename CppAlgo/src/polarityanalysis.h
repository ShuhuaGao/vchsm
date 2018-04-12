#pragma once
#include "PicosElement.h"
#include <vector>
#include "types.h"
#include "StructDefinitions.h"

// function [picos,pol]=polarityanalysis(picos) MATLAB
// Note that the argument picos is passed in and then returned, which implies this function
// will modify this argument. In C++, no need to be so tedious. Just pass the argument by 
// reference and the function can make change to it naturally. 

// pol in MATLAB will be returned
int polarityanalysis(PicosStructArray& picos);


