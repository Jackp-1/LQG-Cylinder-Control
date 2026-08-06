/* This file generated automatically. */ 
/* Do not modify. */ 
#include "udf.h" 
#include "prop.h" 
#include "dpm.h" 
extern DEFINE_EXECUTE_AT_END(Control_loop);
extern DEFINE_PROFILE(update_top, thread, position);
extern DEFINE_PROFILE(update_bot, thread, position);
extern DEFINE_PROFILE(inlet_vel, thread, position);
__declspec(dllexport) UDF_Data udf_data[] = { 
{"Control_loop", (void (*)(void))Control_loop, UDF_TYPE_EXECUTE_AT_END},
{"update_top", (void (*)(void))update_top, UDF_TYPE_PROFILE},
{"update_bot", (void (*)(void))update_bot, UDF_TYPE_PROFILE},
{"inlet_vel", (void (*)(void))inlet_vel, UDF_TYPE_PROFILE},
}; 
__declspec(dllexport) int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data); 
#include "version.h" 
__declspec(dllexport) void UDF_Inquire_Release(int *major, int *minor, int *revision) 
{ 
*major = RampantReleaseMajor; 
*minor = RampantReleaseMinor; 
*revision = RampantReleaseRevision; 
} 
