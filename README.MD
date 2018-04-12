# INTRA-LINGUAL AND CROSS-LINGUAL VOICE CONVERSION USING HARMONIC PLUS STOCHASTIC MODELS
This is a *C++ 11* implementation of the above algorithm proposed by **Dr. Daniel Erro Eslava** based on his PhD thesis, available at [http://www.lsi.upc.edu/~nlp/papers/phd_daniel_erro.pdf]. This algorithm can be used to convert a speaker's voice to another speaker's after training with audio materials from both speakers. For fun, you could convert your voice to President Trump's by downloading his audios from YouTube, for example.

# Features
+ Written in modern C++ (C++11/14)
+ Depend on the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) linear algebra library for high numeric computation performance
+ Utilize OpenMP on Windows or GCD on Mac/iOS to take full advantage of parallel computation on multiple cores of a modern CPU for high speed
+ Expose **C interfaces** to facilitate interoperation with other languages on multiple platforms, such as Windows, Mac OS and iOS

# Usage
Before training and conversion, please first add the proper paths like `/include/vchsm` and `/external/Eigen3.3.4` in your IDE's header file search path such that `train_C.h` and `convert_C.h` can be found and the source files can be compiled correctly.
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
+ For both training and conversion APIs, a client needs to provide multiple arguments such as the list of source audios, target audios and the path of the model file. To make this process more friendly and more concise, the API can also make use of a single configuration file to specify all the required parameters. 
+ **Requirements of input audios:**
	- Mono WAV, i.e., a single channel;
	- 16 bits per sample;
	- 16kHz sampling frequency.
(You may use tools such as [MediaInfo](https://mediaarea.net/en/MediaInfo/Download) to get technical details of an audio file.)
	- For training, the materials for the target speaker and the source speaker should remain the same. For example, in source-1.wav and target-1.wav, they may both say "make America great again". 
	- After training finished, the source speaker can speak anything he/she likes, and the algorithm will convert the voice into one like the target speaker's using the model file obtained in training.

## Example
In the *audios* subdirectory, we have placed the audio samples from *Dr. Daniel Erro Eslava*, which includes a training set composed 20 audio samples from two speakers and a test set containing 10 audio samples from the source speaker.
Let's walk through a simple example. Using the training set in `Audios`, we can train a model with the given 20 audio samples from both the source speaker and target speaker. After that, this generated model can be used to convert any voice record of the source speaker into the voice of the target speaker, i.e., voice conversion.  Therefore, you may make yourself sound like President Trump if you can get his voice samples for training.
On a Intel Core i7 CPU with 4 cores, it takes about 30 seconds to train and a negligible time to convert, which makes this implementation close to *real-time*.

```c++
#include <cstdio>
#include "../include/vchsm/train_C.h"
#include "../include/vchsm/convert_C.h"

#define VERBOSE_TRUE 1
#define VERBOSE_FALSE 0

int main()
{
	// train a model with the given 20 source and target speaker's audios
	// the audios files are named from 1 to 20
	const char* sourceAudioDir = "E:/GitHub/vchsm/Audios/source_train/";
	const char*  targetAudioDir = "E:/GitHub/vchsm/Audios/target_train/";
	const int numTrainSamples = 20;
	const char* sourceAudioList[numTrainSamples];
	const char* targetAudioList[numTrainSamples];
	for (int i = 0; i < numTrainSamples; ++i)
	{
		char* buff = new char[100];
		std::sprintf(buff, "%s%d.wav", sourceAudioDir, i + 1);
		sourceAudioList[i] = buff;
		buff = new char[100];
		std::sprintf(buff, "%s%d.wav", targetAudioDir, i + 1);
		targetAudioList[i] = buff;
	}
	// model file to be generated
	const char* modelFile = "E:/GitHub/vchsm/models/Model.dat";
	// start training	
	trainHSMModel(sourceAudioList, targetAudioList, numTrainSamples, 4, modelFile, VERBOSE_TRUE);
	// deallocate
	for (int i = 0; i < numTrainSamples; ++i)
	{
		delete[] sourceAudioList[i];
		delete[] targetAudioList[i];
	}
	// perform conversion
	const char* testAudio = "E:/GitHub/vchsm/Audios/test/jal_in_42_3.wav";
	const char* testAudioConverted = "E:/GitHub/vchsm/Audios/test/jal_in_42_3_c.wav";
	convertSingle(modelFile, testAudio, testAudioConverted, VERBOSE_TRUE);
	// now we can compare the above audio before and after conversion
	std::getchar();
}
```
Codes for the above example are placed in [./CppAlgo/example/main.cpp](./CppAlgo/example/main.cpp).

## Additional notes for application on Mac/iOS
Since the Mac/iOS platforms prefer GCD to OpenMP for multiple-core parallel programming, the equivalent implementation using GCD is provided in *train_C.mm*.  Just use *train_C.mm*  to replace *train_C.cpp* (which depends on OpenMP) if you plan to deploy it on Mac/iOS. To facilitate the use of this library with Swift, a bridge file *C_Swift_bridge_header.h* is also included. 

# Publication
Wu, Xiaoling, Shuhua Gao, Dong-Yan Huang, and Cheng Xiang. "Voichap: A standalone real-time voice change application on iOS platform." In Asia-Pacific Signal and Information Processing Association Annual Summit and Conference (APSIPA ASC), 2017, pp. 728-732. IEEE, 2017. 
[get pdf here](http://www.apsipa.org/proceedings/2017/CONTENTS/papers2017/14DecThursday/TA-03/TA-03.5.pdf)
# TODO
Add cmake support.