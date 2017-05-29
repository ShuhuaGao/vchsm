#pragma once
#include <Eigen/Dense>
#include <utility>
#include "types.h"

/**
 * Splits signal into overlapped frames using indexing.
 * Implementation of the MATLAB function vec2frames in the source code.
 *@note This function is only called in mfcc.m by "frames = vec2frames( speech, Nw, Ns, 'cols',window, false );". 
 *Therefore, these parameters are fixed: direction = 'cols', window = @hamming and padding = false.
 *Since C++ is not as flexible as MATLAB, for simplicity here we ignore these fixed parameters.
 *Besides, it is found that only the first return value is ever used though two values are returned in vec2frames.m.
*/
Eigen::TMatrixX vec2frames(Eigen::Ref<const Eigen::TVectorX>vec, int Nw, int Ns);