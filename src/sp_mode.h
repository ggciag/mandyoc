#ifndef SP_MODE_H
#define SP_MODE_H

// Modes of surface processes
//
// When adding new modes:
//     1. Add to the enum
//     2. Add to valid_modes[] array
//     2. Add to string_to_sp_mode (in reader.cpp)
//     3. Add to sp_mode_to_string (in reader.cpp)
//

typedef enum {
    SP_NONE,      // Dummy value for inactive state
    SP_DIFFUSION, // "diffusion"
    SP_SEDIMENTATION_ONLY, // "sedimentation_only"
    SP_DIFFUSION_SEDIMENTATION_ONLY , // "diffusion_sedimentation_only"
} SP_Mode;

// Conversion functions
SP_Mode sp_mode_from_string(const char* str);
const char* sp_mode_to_string(SP_Mode mode);

static const char* valid_modes[] = {
    "\"none\"               - No surface processes (when disabled)",
    "\"diffusion\"          - Simple diffusion (requires sp_d_c)",
    "\"sedimentation_only\" - Only sedimentation bellow height adjusted by sea level (requires sea_level)",
    "\"diffusion_sedimentation_only\" - Submarine diffusion-based sedimentation (requires sea_level, sp_Ks, sp_lambda_s)",
    NULL  // Terminator
};

#endif // SP_MODE_H
