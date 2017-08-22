# INTRA-LINGUAL AND CROSS-LINGUAL VOICE CONVERSION USING HARMONIC PLUS STOCHASTIC MODELS
This is a *C++ 11* implementation of the above algorithm proposed by **Dr. Daniel Erro Eslava** based on his PhD thesis, available at [http://www.lsi.upc.edu/~nlp/papers/phd_daniel_erro.pdf]. This algorithm can be used to convert a speaker's voice to another speaker's after training with audio materials from both speakers. For fun, you could convert your voice to President Trump's by downloading his audios from YouTube, for example.

# Features
+ Written in modern C++ (C++11/14)
+ Depend on the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) linear algebra library for high numeric computation performance
+ Utilize OpenMP on Windows or GCD on Mac/iOS to take full advantage of parallel computation on multiple cores of a modern CPU for high speed
+ Expose **C interfaces** to facilitate interoperation with other languages on multiple platforms, such as Windows, Mac OS and iOS

# Usage
## Training
The C API for the training phase are presented in *train_C.h*， which requires the client to specify the source audios, target audios such that it can generate a model file denoting the voice conversion from the source speaker to the target speaker. All the parameters are documented in details in the corresponding header file.
To use an API, just 
```c++
#include "train_C.h"
```
## Conversion
After training, this algorithm produces a model file (whose path is specified by the user in the training phase). Then we can convert any voice audios from the source speaker to make it hear like the target speaker's voice with this model file. All the C APIs for conversion are contained in *convert_C.h,* which are well documented in the header file.
```c++
#include "convert_C.h"
```
## Parallel computing
For the training phase, the algorithm needs to analyze multiple training samples (audios). To speed up the training procedure, a natural way is to divide the training set into several subsets whereby each subset can be processed simultaneously on a multiple-core processor. A naive implementation of this parallel process idea is realized with [OpenMP](https://en.wikipedia.org/wiki/OpenMP) in `train_c.cpp` or [GCD](https://en.wikipedia.org/wiki/Grand_Central_Dispatch) in `train_c.mm`.

## Note
For both training and conversion APIs, a client needs to provide multiple arguments such as the list of source audios, target audios and the path of the model file. To make this process more friendly and more concise, the API can also make use of a single configuration file to specify all the required parameters. 

## Example
In the *audios* subdirectory, we have placed the audio samples from *Dr. Daniel Erro Eslava*, which includes a training set composed 20 audio samples from two speakers and a test set containing 10 audio samples from the source speaker.

## Additional notes for application on Mac/iOS
Since the Mac/iOS platforms prefer GCD to OpenMP for multiple-core parallel programming, the equivalent implementation using GCD is provided in *train_C.mm*.  Just use *train_C.mm*  to replace *train_C.cpp* (which depends on OpenMP) if you plan to deploy it on Mac/iOS. To facilitate the use of this library with Swift, a bridge file *C_Swift_bridge_header.h* is also included. 

