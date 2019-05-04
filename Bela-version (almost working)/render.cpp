#include <Bela.h>
/*
 ____  _____ _        _
 | __ )| ____| |      / \
 |  _ \|  _| | |     / _ \
 | |_) | |___| |___ / ___ \
 |____/|_____|_____/_/   \_\
 
 The platform for ultra-low latency audio and sensor processing
 
 http://bela.io
 
 A project of the Augmented Instruments Laboratory within the
 Centre for Digital Music at Queen Mary University of London.
 http://www.eecs.qmul.ac.uk/~andrewm
 
 (c) 2016 Augmented Instruments Laboratory: Andrew McPherson,
 Astrid Bin, Liam Donovan, Christian Heinrichs, Robert Jack,
 Giulio Moro, Laurel Pardue, Victor Zappi. All rights reserved.
 
 The Bela software is distributed under the GNU Lesser General Public License
 (LGPL 3.0), available here: https://www.gnu.org/licenses/lgpl-3.0.txt
 */

/*
 *    USING A CUSTOM RENDER.CPP FILE FOR PUREDATA PATCHES - LIBPD
 *  ===========================================================
 *  ||                                                       ||
 *  || OPEN THE ENCLOSED _main.pd PATCH FOR MORE INFORMATION ||
 *  || ----------------------------------------------------- ||
 *  ===========================================================
 */
#include <Bela.h>
#include <DigitalChannelManager.h>
#include <cmath>
#include <stdio.h>
#define PD_THREADED_IO
#include <libpd/z_libpd.h>
extern "C" {
#include <libpd/s_stuff.h>
};
#include <UdpServer.h>
#include <Midi.h>
#include <Scope.h>
#include <string>
#include <sstream>
#include <algorithm>
// my libs
#include <math.h>
#include <vector>
#include <iostream>

enum { minFirstDigitalChannel = 10 };
static unsigned int gAnalogChannelsInUse;
static unsigned int gDigitalChannelsInUse;
static unsigned int gScopeChannelsInUse = 4;
static unsigned int gLibpdBlockSize;
static unsigned int gChannelsInUse;
//static const unsigned int gFirstAudioChannel = 0;
static unsigned int gFirstAnalogInChannel;
static unsigned int gFirstAnalogOutChannel;
static unsigned int gFirstDigitalChannel;
static unsigned int gLibpdDigitalChannelOffset;
static unsigned int gFirstScopeChannel;

void Bela_userSettings(BelaInitSettings *settings)
{
    settings->uniformSampleRate = 1;
    settings->interleave = 0;
    settings->analogOutputsPersist = 0;
}

/*
 *  MODIFICATION
 *  ------------
 *  Global variables for tremolo effect applied to libpd output
 */

//float gTremoloRate = 4.0;
//float gPhase;
//float gInverseSampleRate;

// function declaration
double fireflySimulation(int f, double samplerate);
double median(std::vector<double> medi);

// global intialisation
//int in = 0;
//double interactions [1][11]; //interactions[0] = {1,2,3,4,5,6,7,8,9,10,11};
std::vector<double> E;
std::vector<double> C;
std::vector<double> domega(2,0);
int fire = 0;

int simlength = 8e4;
int sr =1e3;

// initial phase (temp2 = rand(2,1);)
std::vector<double> temp2={0.191519450378892, 0.622108771039832};
double pas_phi = temp2[0];
std::vector<double> phi_ext;

// initial frequency (temp = 2*2.^((rand(2,1)-0.5)*2);)
std::vector<double> temp={1.834587198568937,2.970523429156172};
double omega = temp[0];
//omega_ext= temp[1];
int f_ = 0;
//int ii;
//bool count = true;


#define pi 3.14159265358979323846


/*********/

float* gInBuf;
float* gOutBuf;
#define PARSE_MIDI
static std::vector<Midi*> midi;
std::vector<std::string> gMidiPortNames;

void dumpMidi()
{
    if(midi.size() == 0)
    {
        printf("No MIDI device enabled\n");
        return;
    }
    printf("The following MIDI devices are enabled:\n");
    printf("%4s%20s %3s %3s %s\n",
           "Num",
           "Name",
           "In",
           "Out",
           "Pd channels"
           );
    for(unsigned int n = 0; n < midi.size(); ++n)
    {
        printf("[%2d]%20s %3s %3s (%d-%d)\n",
               n,
               gMidiPortNames[n].c_str(),
               midi[n]->isInputEnabled() ? "x" : "_",
               midi[n]->isOutputEnabled() ? "x" : "_",
               n * 16 + 1,
               n * 16 + 16
               );
    }
}

Midi* openMidiDevice(std::string name, bool verboseSuccess = false, bool verboseError = false)
{
    Midi* newMidi;
    newMidi = new Midi();
    newMidi->readFrom(name.c_str());
    newMidi->writeTo(name.c_str());
#ifdef PARSE_MIDI
    newMidi->enableParser(true);
#else
    newMidi->enableParser(false);
#endif /* PARSE_MIDI */
    if(newMidi->isOutputEnabled())
    {
        if(verboseSuccess)
            printf("Opened MIDI device %s as output\n", name.c_str());
    }
    if(newMidi->isInputEnabled())
    {
        if(verboseSuccess)
            printf("Opened MIDI device %s as input\n", name.c_str());
    }
    if(!newMidi->isInputEnabled() && !newMidi->isOutputEnabled())
    {
        if(verboseError)
            fprintf(stderr, "Failed to open  MIDI device %s\n", name.c_str());
        return nullptr;
    } else {
        return newMidi;
    }
}

static unsigned int getPortChannel(int* channel){
    unsigned int port = 0;
    while(*channel > 16){
        *channel -= 16;
        port += 1;
    }
    if(port >= midi.size()){
        // if the port number exceeds the number of ports available, send out
        // of the first port
        rt_fprintf(stderr, "Port out of range, using port 0 instead\n");
        port = 0;
    }
    return port;
}

void Bela_MidiOutNoteOn(int channel, int pitch, int velocity) {
    int port = getPortChannel(&channel);
    rt_printf("noteout _ port: %d, channel: %d, pitch: %d, velocity %d\n", port, channel, pitch, velocity);
    midi[port]->writeNoteOn(channel, pitch, velocity);
}

void Bela_MidiOutControlChange(int channel, int controller, int value) {
    int port = getPortChannel(&channel);
    rt_printf("ctlout _ port: %d, channel: %d, controller: %d, value: %d\n", port, channel, controller, value);
    midi[port]->writeControlChange(channel, controller, value);
}

void Bela_MidiOutProgramChange(int channel, int program) {
    int port = getPortChannel(&channel);
    rt_printf("pgmout _ port: %d, channel: %d, program: %d\n", port, channel, program);
    midi[port]->writeProgramChange(channel, program);
}

void Bela_MidiOutPitchBend(int channel, int value) {
    int port = getPortChannel(&channel);
    rt_printf("bendout _ port: %d, channel: %d, value: %d\n", port, channel, value);
    midi[port]->writePitchBend(channel, value);
}

void Bela_MidiOutAftertouch(int channel, int pressure){
    int port = getPortChannel(&channel);
    rt_printf("touchout _ port: %d, channel: %d, pressure: %d\n", port, channel, pressure);
    midi[port]->writeChannelPressure(channel, pressure);
}

void Bela_MidiOutPolyAftertouch(int channel, int pitch, int pressure){
    int port = getPortChannel(&channel);
    rt_printf("polytouchout _ port: %d, channel: %d, pitch: %d, pressure: %d\n", port, channel, pitch, pressure);
    midi[port]->writePolyphonicKeyPressure(channel, pitch, pressure);
}

void Bela_MidiOutByte(int port, int byte){
    rt_printf("port: %d, byte: %d\n", port, byte);
    if(port > (int)midi.size()){
        // if the port is out of range, redirect to the first port.
        rt_fprintf(stderr, "Port out of range, using port 0 instead\n");
        port = 0;
    }
    midi[port]->writeOutput(byte);
}

void Bela_printHook(const char *received){
    rt_printf("%s", received);
}

static DigitalChannelManager dcm;

void sendDigitalMessage(bool state, unsigned int delay, void* receiverName){
    libpd_float((const char*)receiverName, (float)state);
    //    rt_printf("%s: %d\n", (char*)receiverName, state);
}

void Bela_messageHook(const char *source, const char *symbol, int argc, t_atom *argv){
    if(strcmp(source, "bela_setMidi") == 0){
        int num[3] = {0, 0, 0};
        for(int n = 0; n < argc && n < 3; ++n)
        {
            if(!libpd_is_float(&argv[n]))
            {
                fprintf(stderr, "Wrong format for Bela_setMidi, expected:[hw 1 0 0(");
                return;
            }
            num[n] = libpd_get_float(&argv[n]);
        }
        std::ostringstream deviceName;
        deviceName << symbol << ":" << num[0] << "," << num[1] << "," << num[2];
        printf("Adding Midi device: %s\n", deviceName.str().c_str());
        Midi* newMidi = openMidiDevice(deviceName.str(), false, true);
        if(newMidi)
        {
            midi.push_back(newMidi);
            gMidiPortNames.push_back(deviceName.str());
        }
        dumpMidi();
        return;
    }
    if(strcmp(source, "bela_setDigital") == 0){
        // symbol is the direction, argv[0] is the channel, argv[1] (optional)
        // is signal("sig" or "~") or message("message", default) rate
        bool isMessageRate = true; // defaults to message rate
        bool direction = 0; // initialize it just to avoid the compiler's warning
        bool disable = false;
        if(strcmp(symbol, "in") == 0){
            direction = INPUT;
        } else if(strcmp(symbol, "out") == 0){
            direction = OUTPUT;
        } else if(strcmp(symbol, "disable") == 0){
            disable = true;
        } else {
            return;
        }
        if(argc == 0){
            return;
        } else if (libpd_is_float(&argv[0]) == false){
            return;
        }
        int channel = libpd_get_float(&argv[0]) - gLibpdDigitalChannelOffset;
        if(disable == true){
            dcm.unmanage(channel);
            return;
        }
        if(argc >= 2){
            t_atom* a = &argv[1];
            if(libpd_is_symbol(a)){
                char *s = libpd_get_symbol(a);
                if(strcmp(s, "~") == 0  || strncmp(s, "sig", 3) == 0){
                    isMessageRate = false;
                }
            }
        }
        dcm.manage(channel, direction, isMessageRate);
        return;
    }
}

void Bela_floatHook(const char *source, float value){
    /*
     *  MODIFICATION
     *  ------------
     *  Parse float sent to receiver 'tremoloRate' and assign it to a global variable
     *  N.B. When using libpd receiver names need to be registered (see setup() function below)
     */
    //    if(strncmp(source, "tremoloRate", 11) == 0){
    //        gTremoloRate = value;
    //    }
    
    /*********/
    
    // let's make this as optimized as possible for built-in digital Out parsing
    // the built-in digital receivers are of the form "bela_digitalOutXX" where XX is between gLibpdDigitalChannelOffset and (gLibpdDigitalCHannelOffset+gDigitalChannelsInUse)
    static int prefixLength = 15; // strlen("bela_digitalOut")
    if(strncmp(source, "bela_digitalOut", prefixLength)==0){
        if(source[prefixLength] != 0){ //the two ifs are used instead of if(strlen(source) >= prefixLength+2)
            if(source[prefixLength + 1] != 0){
                // quickly convert the suffix to integer, assuming they are numbers, avoiding to call atoi
                int receiver = ((source[prefixLength] - 48) * 10);
                receiver += (source[prefixLength+1] - 48);
                unsigned int channel = receiver - gLibpdDigitalChannelOffset; // go back to the actual Bela digital channel number
                if(channel < gDigitalChannelsInUse){ //number of digital channels
                    dcm.setValue(channel, value);
                }
            }
        }
    }
}



std::vector<std::string> gReceiverInputNames;
std::vector<std::string> gReceiverOutputNames;
void generateDigitalNames(unsigned int numDigitals, unsigned int libpdOffset, std::vector<std::string>& receiverInputNames, std::vector<std::string>& receiverOutputNames)
{
    std::string inBaseString = "bela_digitalIn";
    std::string outBaseString = "bela_digitalOut";
    for(unsigned int i = 0; i<numDigitals; i++)
    {
        receiverInputNames.push_back(inBaseString + std::to_string(i+libpdOffset));
        receiverOutputNames.push_back(outBaseString + std::to_string(i+libpdOffset));
    }
}

void printDigitalNames(std::vector<std::string>& receiverInputNames, std::vector<std::string>& receiverOutputNames)
{
    printf("DIGITAL INPUTS\n");
    for(unsigned int i=0; i<gDigitalChannelsInUse; i++)
        printf("%s\n", receiverInputNames[i].c_str());
    printf("DIGITAL OUTPUTS\n");
    for(unsigned int i=0; i<gDigitalChannelsInUse; i++)
        printf("%s\n", receiverOutputNames[i].c_str());
}

static char multiplexerArray[] = {"bela_multiplexer"};
static int multiplexerArraySize = 0;
static bool pdMultiplexerActive = false;

#ifdef PD_THREADED_IO
void fdLoop(void* arg){
    // t_pdinstance* pd_that = (t_pdinstance*)arg;
    while(!gShouldStop){
        sys_doio();
        //sys_doio(pd_that);
        usleep(3000);
    }
}
#endif /* PD_THREADED_IO */

Scope scope;
float* gScopeOut;
void* gPatch;
bool gDigitalEnabled = 0;

bool setup(BelaContext *context, void *userData)
{
    // Check Pd's version
    int major, minor, bugfix;
    sys_getversion(&major, &minor, &bugfix);
    printf("Running Pd %d.%d-%d\n", major, minor, bugfix);
    // We requested in Bela_userSettings() to have uniform sampling rate for audio
    // and analog and non-interleaved buffers.
    // So let's check this actually happened
    if(context->analogSampleRate != context->audioSampleRate)
    {
        fprintf(stderr, "The sample rate of analog and audio must match. Try running with --uniform-sample-rate\n");
        return false;
    }
    if(context->flags & BELA_FLAG_INTERLEAVED)
    {
        fprintf(stderr, "The audio and analog channels must be interleaved.\n");
        return false;
    }
    
    if(context->digitalFrames > 0 && context->digitalChannels > 0)
        gDigitalEnabled = 1;
    
    /*
     *  MODIFICATION
     *  ------------
     *  Initialise variables for tremolo effect
     */
    
    //    gInverseSampleRate = 1.0 / context->audioSampleRate;
    //    gPhase = 0.0;
    
    /*********/
    
    // add here other devices you need
    gMidiPortNames.push_back("hw:1,0,0");
    //gMidiPortNames.push_back("hw:0,0,0");
    //gMidiPortNames.push_back("hw:1,0,1");
    
    scope.setup(gScopeChannelsInUse, context->audioSampleRate);
    gScopeOut = new float[gScopeChannelsInUse];
    
    // Check first of all if the patch file exists. Will actually open it later.
    char file[] = "_main.pd";
    char folder[] = "./";
    unsigned int strSize = strlen(file) + strlen(folder) + 1;
    char* str = (char*)malloc(sizeof(char) * strSize);
    snprintf(str, strSize, "%s%s", folder, file);
    if(access(str, F_OK) == -1 ) {
        printf("Error file %s/%s not found. The %s file should be your main patch.\n", folder, file, file);
        return false;
    }
    free(str);
    
    // analog setup
    gAnalogChannelsInUse = context->analogInChannels;
    gDigitalChannelsInUse = context->digitalChannels;
    printf("Audio channels in use: %d\n", context->audioOutChannels);
    printf("Analog channels in use: %d\n", gAnalogChannelsInUse);
    printf("Digital channels in use: %d\n", gDigitalChannelsInUse);
    
    // Channel distribution
    gFirstAnalogInChannel = std::max(context->audioInChannels, context->audioOutChannels);
    gFirstAnalogOutChannel = gFirstAnalogInChannel;
    gFirstDigitalChannel = gFirstAnalogInChannel + std::max(context->analogInChannels, context->analogOutChannels);
    if(gFirstDigitalChannel < minFirstDigitalChannel)
        gFirstDigitalChannel = minFirstDigitalChannel; //for backwards compatibility
    gLibpdDigitalChannelOffset = gFirstDigitalChannel + 1;
    gFirstScopeChannel = gFirstDigitalChannel + gDigitalChannelsInUse;
    
    gChannelsInUse = gFirstScopeChannel + gScopeChannelsInUse;
    
    // Create receiverNames for digital channels
    generateDigitalNames(gDigitalChannelsInUse, gLibpdDigitalChannelOffset, gReceiverInputNames, gReceiverOutputNames);
    
    // digital setup
    if(gDigitalEnabled)
    {
        dcm.setCallback(sendDigitalMessage);
        if(gDigitalChannelsInUse > 0){
            for(unsigned int ch = 0; ch < gDigitalChannelsInUse; ++ch){
                dcm.setCallbackArgument(ch, (void*) gReceiverInputNames[ch].c_str());
            }
        }
    }
    
    unsigned int n = 0;
    while(n < gMidiPortNames.size())
    {
        Midi* newMidi = openMidiDevice(gMidiPortNames[n], false, false);
        if(newMidi)
        {
            midi.push_back(newMidi);
            ++n;
        } else {
            gMidiPortNames.erase(gMidiPortNames.begin() + n);
        }
    }
    dumpMidi();
    
    // check that we are not running with a blocksize smaller than gLibPdBlockSize
    gLibpdBlockSize = libpd_blocksize();
    if(context->audioFrames < gLibpdBlockSize){
        fprintf(stderr, "Error: minimum block size must be %d\n", gLibpdBlockSize);
        return false;
    }
    
    // set hooks before calling libpd_init
    libpd_set_printhook(Bela_printHook);
    libpd_set_floathook(Bela_floatHook);
    libpd_set_messagehook(Bela_messageHook);
    libpd_set_noteonhook(Bela_MidiOutNoteOn);
    libpd_set_controlchangehook(Bela_MidiOutControlChange);
    libpd_set_programchangehook(Bela_MidiOutProgramChange);
    libpd_set_pitchbendhook(Bela_MidiOutPitchBend);
    libpd_set_aftertouchhook(Bela_MidiOutAftertouch);
    libpd_set_polyaftertouchhook(Bela_MidiOutPolyAftertouch);
    libpd_set_midibytehook(Bela_MidiOutByte);
    
    //initialize libpd. This clears the search path
    libpd_init();
    //Add the current folder to the search path for externals
    libpd_add_to_search_path(".");
    libpd_add_to_search_path("../pd-externals");
    
    libpd_init_audio(gChannelsInUse, gChannelsInUse, context->audioSampleRate);
    gInBuf = get_sys_soundin();
    gOutBuf = get_sys_soundout();
    
    // start DSP:
    // [; pd dsp 1(
    libpd_start_message(1);
    libpd_add_float(1.0f);
    libpd_finish_message("pd", "dsp");
    
    // Bind your receivers here
    for(unsigned int i = 0; i < gDigitalChannelsInUse; i++)
        libpd_bind(gReceiverOutputNames[i].c_str());
    libpd_bind("bela_setDigital");
    libpd_bind("bela_setMidi");
    
    /*
     *  MODIFICATION
     *  ------------
     *  Bind an additional receiver for the tremoloRate parameter
     */
    //    libpd_bind("tremoloRate");
    
    
    // if length of phase > 3333 -> phase(end) = phase(0)
    
    //    if (E.size()==0){
    //    local intialisation
    double omega_ext = temp[1];
    //        int f = 0;
    std::vector<double> omegas;
    std::vector<double> phase;
    
    phi_ext.push_back(temp2[1]);
    E.push_back(1);
    C.push_back(1);
    
    phase.push_back(pas_phi);
    omegas.push_back(omega);
    
    //        count==false;
    //    }
    
    //for (int i=2-1; i<simlength; i++) {
    //        ii=i;
    if (phi_ext[phi_ext.size()-1]<1) {
        phi_ext.push_back(phi_ext[phi_ext.size()-1]+omega_ext/sr); //tooth saw -> increase (phase=amplitude)
        if (phi_ext[phi_ext.size()-1]>=1){
            phi_ext[phi_ext.size()-1]=1;
        }
    }
    else{
        f_=(f_+1)%2;
        // std::cout << f_;
        phi_ext.push_back(0);
    }
    phase.push_back(fireflySimulation(f_,sr));
    //    omegas.push_back(omega);
    libpd_start_message(1);
    libpd_add_float(f_);
    libpd_finish_list("firefly");
    //}
    
    libpd_bang("foo");
    
    libpd_start_message(1);
    libpd_add_float(1);
    libpd_finish_list("test");
    
    
    
    /*********/
    
    // open patch:
    gPatch = libpd_openfile(file, folder);
    if(gPatch == NULL){
        printf("Error: file %s/%s is corrupted.\n", folder, file);
        return false;
    }
    
    // If the user wants to use the multiplexer capelet,
    // the patch will have to contain an array called "bela_multiplexer"
    // and a receiver [r bela_multiplexerChannels]
    if(context->multiplexerChannels > 0 && libpd_arraysize(multiplexerArray) >= 0){
        pdMultiplexerActive = true;
        multiplexerArraySize = context->multiplexerChannels * context->analogInChannels;
        // [; bela_multiplexer ` multiplexerArraySize` resize(
        libpd_start_message(1);
        libpd_add_float(multiplexerArraySize);
        libpd_finish_message(multiplexerArray, "resize");
        // [; bela_multiplexerChannels `context->multiplexerChannels`(
        libpd_float("bela_multiplexerChannels", context->multiplexerChannels);
    }
    
    // Tell Pd that we will manage the io loop,
    // and we do so in an Auxiliary Task
#ifdef PD_THREADED_IO
    sys_dontmanageio(1);
    AuxiliaryTask fdTask;
    fdTask = Bela_createAuxiliaryTask(fdLoop, 50, "libpd-fdTask", NULL);
    Bela_scheduleAuxiliaryTask(fdTask);
#endif /* PD_THREADED_IO */
    
    //    dcm.setVerbose(false);
    return true;
}



double fireflySimulation(int f, double sr) {
    // if length of phase > 3333 -> phase(end) = phase(0)
    
    // local variable declaration
    const double Alpha = 0.4;
    const double Beta = 0.4;
    const int FilterLength = 8;
    const double tref = 5;
    const int DoublingThreshold = 20;
    double phi;
    double phiref;
    std::vector<double> CE;
    
    if (pas_phi < 1){
        phi = pas_phi+omega/sr; //tooth saw -> increase (phase=amplitude)
        if (phi>=1) {phi=1;}
    }
    else{
        fire = (fire+1)%2; //fire only on every other peak
        phi = 0;
    }
    
    phiref = tref*omega*5e-4;
    
    if (phi_ext[phi_ext.size()-1] >= 1){    // if external node is at maximum
        if (f == 0){
            //            pas_.in=pas_.in+1;
            //            pas_.interactions(pas_.in,1:11) = [-1 i j 0 0 0 0 0 0 0 0];
        }
        if (f == 1 && fire == 1){ // node j pas_.fires and firing has already occured at this time.
            //        pas_.in=pas_.in+1;
            //        pas_.interactions(pas_.in,1:3) = [0 i j];
        }
        
        if (f == 1 && fire == 0){// if node j pas_.fires and firing has not already occured at this time.
            //            phi_ext[phi_ext.size()-1]  = 1;
            
            if (phi > phiref){ // if the phase of node k is above phiref
                
                // PHASE COUPLING FUNCTION
                phi = phi+Alpha*(-sin(phi*2*pi))*std::abs(sin(phi*2*pi));
                if (phi<0){phi=0;} // trim phase to [0 1] range
                //                else if (phi>=1){phi=1;}
                if (phi>=1){phi=1;}
                
                if (phi < phiref || phi > 1-phiref){
                    E.push_back(0);
                }
                else{
                    E.push_back((1-cos(2*pi*phi))/2); //else, error is calculated like this
                }
                
                if (E.size() < FilterLength){
                    C.push_back(median(E));
                }
                else{
                    for(int ii =0;ii<(FilterLength);ii++){
                        CE.push_back(E[E.size()-FilterLength+ii]);
                    }
                    C.push_back(median(CE));
                }
                
                // FREQUENCY COUPLING FUNCTION
                domega[0] = ((omega * pow(2,Beta*-sin(2*pi*phi)*C[C.size()-1]))+domega[0]*domega[1])/(domega[1]+1);
                domega[1] = domega[1]+1;
                if (domega[1] > DoublingThreshold){
                    omega = 2*omega;
                    domega = {0,0};
                }
                //            pas_.in=pas_.in+1;
                //            pas_.interactions(pas_.in,1:9) = [1 i j k 1 pas_.C(end) phi-pas_.phi Beta*sin(2*pi*phi) pas_.domega(1,1)];
            }
            else{
                //            pas_.in=pas_.in+1;
                //            pas_.interactions(pas_.in,1:6) = [1 i j k 0 pas_.C(end)]; %col 5 == 0 means pas_.fire within tref
            }
        }
        
        if (domega[1] > 0){ // if there has been freq adjustments
            omega = domega[0];    // Reachback pas_.firefly algorithm
            domega = {0,0};
            //                   if pas_.omega(j) == 0
            //                                 disp('a node has died... pas_.omega = zero')
            //                         end
            //        pas_.in=pas_.in+1;
            //        pas_.interactions(pas_.in,1:10) = [-2 i j j 1 pas_.C(end) -1 0 0 pas_.omega(j)];
        }
        
        
        //    if (length(pas_.E{j}) < FilterLength
        //        allzeroes = sum(pas_.E{j}(1:length(pas_.E{j})));
        //    allzeroes = filterType(pas_.E{j}(1:length(pas_.E{j})));
        //    else{
        //        allzeroes = sum(pas_.E{j}(length(pas_.E{j})-FilterLength+1:length(pas_.E{j})));
        //    allzeroes = filterType(pas_.E{j}(length(pas_.E{j})-FilterLength+1:length(pas_.E{j})));
        //        }
        //                 if allzeroes == 0 && length(pas_.interactions(pas_.interactions(:,3)==j & pas_.interactions(:,11)==1,1)) > 8
        //                        % length(pas_.interactions(pas_.interactions(:,3)==j,1))
        //                        % disp('selfmean adjust')
        //                        sorted.Omegas = sort(pas_.interactions(pas_.interactions(:,3)==j & pas_.interactions(:,11)==1,2),'descend')';
        //
        //                        pas_.omega(j) = -2000/mean(diff(sorted.Omegas(1:FilterLength)));
        //                    pas_.in=pas_.in+1;
        //                    pas_.interactions(pas_.in,1:10) = [-3 i j j 1 pas_.C(end) -1 0 0 pas_.omega(j)];
        //                end
        //    end
        
    }
    if (phi>1){phi=1;}
    
    pas_phi = phi;
    
    return phi;
}




double median(std::vector<double> medi)
{
    std::sort(medi.begin(), medi.end());     // sort values
    
    double tmedian;
    if (medi.size() % 2 == 0)           // even
        tmedian = (medi[medi.size() / 2 - 1] + medi[medi.size() / 2]) / 2;
    else                                // odd
        tmedian = medi[medi.size() / 2];
    
    return tmedian;
}




void render(BelaContext *context, void *userData)
{
    int num;
    
#ifdef PARSE_MIDI
    for(unsigned int port = 0; port < midi.size(); ++port){
        while((num = midi[port]->getParser()->numAvailableMessages()) > 0){
            static MidiChannelMessage message;
            message = midi[port]->getParser()->getNextChannelMessage();
            rt_printf("On port %d (%s): ", port, gMidiPortNames[port].c_str());
            message.prettyPrint(); // use this to print beautified message (channel, data bytes)
            switch(message.getType()){
                case kmmNoteOn:
                {
                    int noteNumber = message.getDataByte(0);
                    int velocity = message.getDataByte(1);
                    int channel = message.getChannel();
                    libpd_noteon(channel + port * 16, noteNumber, velocity);
                    break;
                }
                case kmmNoteOff:
                {
                    /* PureData does not seem to handle noteoff messages as per the MIDI specs,
                     * so that the noteoff velocity is ignored. Here we convert them to noteon
                     * with a velocity of 0.
                     */
                    int noteNumber = message.getDataByte(0);
                    //                int velocity = message.getDataByte(1); // would be ignored by Pd
                    int channel = message.getChannel();
                    libpd_noteon(channel + port * 16, noteNumber, 0);
                    break;
                }
                case kmmControlChange:
                {
                    int channel = message.getChannel();
                    int controller = message.getDataByte(0);
                    int value = message.getDataByte(1);
                    libpd_controlchange(channel + port * 16, controller, value);
                    break;
                }
                case kmmProgramChange:
                {
                    int channel = message.getChannel();
                    int program = message.getDataByte(0);
                    libpd_programchange(channel + port * 16, program);
                    break;
                }
                case kmmPolyphonicKeyPressure:
                {
                    int channel = message.getChannel();
                    int pitch = message.getDataByte(0);
                    int value = message.getDataByte(1);
                    libpd_polyaftertouch(channel + port * 16, pitch, value);
                    break;
                }
                case kmmChannelPressure:
                {
                    int channel = message.getChannel();
                    int value = message.getDataByte(0);
                    libpd_aftertouch(channel + port * 16, value);
                    break;
                }
                case kmmPitchBend:
                {
                    int channel = message.getChannel();
                    int value =  ((message.getDataByte(1) << 7)| message.getDataByte(0)) - 8192;
                    libpd_pitchbend(channel + port * 16, value);
                    break;
                }
                case kmmSystem:
                    // currently Bela only handles sysrealtime, and it does so pretending it is a channel message with no data bytes, so we have to re-assemble the status byte
                {
                    int channel = message.getChannel();
                    int status = message.getStatusByte();
                    int byte = channel | status;
                    libpd_sysrealtime(port, byte);
                    break;
                }
                case kmmNone:
                case kmmAny:
                    break;
            }
        }
    }
#else
    int input;
    for(unsigned int port = 0; port < NUM_MIDI_PORTS; ++port){
        while((input = midi[port].getInput()) >= 0){
            libpd_midibyte(port, input);
        }
    }
#endif /* PARSE_MIDI */
    unsigned int numberOfPdBlocksToProcess = context->audioFrames / gLibpdBlockSize;
    
    // Remember: we have non-interleaved buffers and the same sampling rate for
    // analogs, audio and digitals
    for(unsigned int tick = 0; tick < numberOfPdBlocksToProcess; ++tick)
    {
        //audio input
        for(int n = 0; n < context->audioInChannels; ++n)
        {
            memcpy(
                   gInBuf + n * gLibpdBlockSize,
                   context->audioIn + tick * gLibpdBlockSize + n * context->audioFrames,
                   sizeof(context->audioIn[0]) * gLibpdBlockSize
                   );
        }
        
        // analog input
        for(int n = 0; n < context->analogInChannels; ++n)
        {
            memcpy(
                   gInBuf + gLibpdBlockSize * gFirstAnalogInChannel + n * gLibpdBlockSize,
                   context->analogIn + tick * gLibpdBlockSize + n * context->analogFrames,
                   sizeof(context->analogIn[0]) * gLibpdBlockSize
                   );
        }
        // multiplexed analog input
        if(pdMultiplexerActive)
        {
            // we do not disable regular analog inputs if muxer is active, because user may have bridged them on the board and
            // they may be using half of them at a high sampling-rate
            static int lastMuxerUpdate = 0;
            if(++lastMuxerUpdate == multiplexerArraySize){
                lastMuxerUpdate = 0;
                libpd_write_array(multiplexerArray, 0, (float *const)context->multiplexerAnalogIn, multiplexerArraySize);
            }
        }
        
        unsigned int digitalFrameBase = gLibpdBlockSize * tick;
        unsigned int j;
        unsigned int k;
        float* p0;
        float* p1;
        // digital input
        if(gDigitalEnabled)
        {
            // digital in at message-rate
            dcm.processInput(&context->digital[digitalFrameBase], gLibpdBlockSize);
            
            // digital in at signal-rate
            for (j = 0, p0 = gInBuf; j < gLibpdBlockSize; j++, p0++) {
                unsigned int digitalFrame = digitalFrameBase + j;
                for (k = 0, p1 = p0 + gLibpdBlockSize * gFirstDigitalChannel;
                     k < 16; ++k, p1 += gLibpdBlockSize) {
                    if(dcm.isSignalRate(k) && dcm.isInput(k)){ // only process input channels that are handled at signal rate
                        *p1 = digitalRead(context, digitalFrame, k);
                    }
                }
            }
        }
        
        libpd_process_sys(); // process the block
        
        // digital outputs
        if(gDigitalEnabled)
        {
            // digital out at signal-rate
            for (j = 0, p0 = gOutBuf; j < gLibpdBlockSize; ++j, ++p0) {
                unsigned int digitalFrame = (digitalFrameBase + j);
                for (k = 0, p1 = p0  + gLibpdBlockSize * gFirstDigitalChannel;
                     k < context->digitalChannels; k++, p1 += gLibpdBlockSize)
                {
                    if(dcm.isSignalRate(k) && dcm.isOutput(k)){ // only process output channels that are handled at signal rate
                        digitalWriteOnce(context, digitalFrame, k, *p1 > 0.5);
                    }
                }
            }
            
            // digital out at message-rate
            dcm.processOutput(&context->digital[digitalFrameBase], gLibpdBlockSize);
        }
        
        // scope output
        for (j = 0, p0 = gOutBuf; j < gLibpdBlockSize; ++j, ++p0) {
            for (k = 0, p1 = p0 + gLibpdBlockSize * gFirstScopeChannel; k < gScopeChannelsInUse; k++, p1 += gLibpdBlockSize) {
                gScopeOut[k] = *p1;
            }
            scope.log(gScopeOut[0], gScopeOut[1], gScopeOut[2], gScopeOut[3]);
        }
        
        /*
         *  MODIFICATION
         *  ------------
         *  Processing for tremolo effect while writing libpd
         *  output to Bela output buffer
         */
        
        //for (j = 0, p0 = gOutBuf; j < gLibpdBlockSize; j++, p0++) {
        //    // Generate a sinewave with frequency set by gTremoloRate
        //    // and amplitude from -0.5 to 0.5
        //    float lfo = sinf(gPhase) * 0.5;
        //    // Keep track and wrap the phase of the sinewave
        //    gPhase += 2.0 * M_PI * gTremoloRate * gInverseSampleRate;
        //    if(gPhase > 2.0 * M_PI)
        //        gPhase -= 2.0 * M_PI;
        //    for (k = 0, p1 = p0; k < context->audioOutChannels; k++, p1 += gLibpdBlockSize) {
        //        // *p1 here is the sample in Pd's buffer which
        //        // corresponds to the jth frame of the kth channel
        //        // we edit its value in place
        //        *p1 = *p1 * lfo;
        //    }
        //}
        
        /*********/
        
        
        //    if (E.size()==0){
        //    local intialisation
        double omega_ext = temp[1];
        //        int f = 0;
        std::vector<double> omegas;
        std::vector<double> phase;
        
        phi_ext.push_back(temp2[1]);
        E.push_back(1);
        C.push_back(1);
        
        phase.push_back(pas_phi);
        omegas.push_back(omega);
        int sou_=0;
        
        //        count==false;
        //    }
        
        //for (int i=2-1; i<simlength; i++) {
        //        ii=i;
        if (phi_ext[phi_ext.size()-1]<1) {
            phi_ext.push_back(phi_ext[phi_ext.size()-1]+omega_ext/sr); //tooth saw -> increase (phase=amplitude)
            if (phi_ext[phi_ext.size()-1]>=1){
                phi_ext[phi_ext.size()-1]=1;
                sou_=1;
            }
        }
        else{
            f_=(f_+1)%2;
            // std::cout << f_;
            phi_ext.push_back(0);
            sou_=0;
        }
        phase.push_back(fireflySimulation(f_,sr));
        //    omegas.push_back(omega);
        
        libpd_start_message(1);
        libpd_add_float(phi_ext[phi_ext.size()-1]);
        //libpd_add_float(sou_);
        libpd_finish_list("firefly");
        //}
        
        //libpd_bang("foo");
        
        //libpd_start_message(1);
        //libpd_add_float(1);
        //libpd_finish_list("test");
        
        
        // audio output
        for(int n = 0; n < context->audioOutChannels; ++n)
        {
            memcpy(
                   context->audioOut + tick * gLibpdBlockSize + n * context->audioFrames,
                   gOutBuf + n * gLibpdBlockSize,
                   sizeof(context->audioOut[0]) * gLibpdBlockSize
                   );
        }
        
        //analog output
        for(int n = 0; n < context->analogOutChannels; ++n)
        {
            memcpy(
                   context->analogOut + tick * gLibpdBlockSize + n * context->analogFrames,
                   gOutBuf + gLibpdBlockSize * gFirstAnalogOutChannel + n * gLibpdBlockSize,
                   sizeof(context->analogOut[0]) * gLibpdBlockSize
                   );
        }
    }
}

void cleanup(BelaContext *context, void *userData)
{
    for(auto a : midi)
    {
        delete a;
    }
    libpd_closefile(gPatch);
    delete [] gScopeOut;
}
