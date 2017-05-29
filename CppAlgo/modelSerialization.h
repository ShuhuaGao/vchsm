#pragma once
#include <string>
#include "StructDefinitions.h"

/**
 * Serialize the model into the file
 *@param model the HSM model
 *@param file a binary file for serialization
*/
void serializeModel(const HSMModel& model, const std::string& file);


/**
 * Recover the HSM model from a binary file.
 *@param file the file containing the model parameters
 *@return the reconstructed HSM model
*/
HSMModel deserializeModel(const std::string& file);

