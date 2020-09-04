varDict_DriftModel = {
    '^[Cc]oncentration': 'concentration',
    '^[Pp]ressure': 'pressure',
    '^[Vv]elocity_*[Xx]': 'velocity_x',
    '^[Vv]elocity_*[Yy]': 'velocity_y',
    '^[Vv]elocity_*[zZ]': 'velocity_z',
    '^[Vv]oid_*[Ff]raction': 'void_fraction',
    '^[Tt]emperature': 'temperature',
    '^[sS]team_*[Ee]nthalpy': 'steam_enthalpy',
    '^[lL]iquid_*[eE]nthalpy': 'liquid_enthalpy',
    '^[Ee]nthalpy': 'enthalpy',
    '^[sS]team_*[Dd]ensity': 'steam_density',
    '^[lL]iquid_*[Dd]ensity': 'liquid_density',
    '^[d]Density': 'density',
}

varDict_SinglePhase = {
    '^[Pp]ressure': 'pressure',
    '^[Vv]elocity_*[Xx]': 'velocity_x',
    '^[Vv]elocity_*[Yy]': 'velocity_y',
    '^[Vv]elocity_*[zZ]': 'velocity_z',
    '^[Tt]emperature': 'temperature',
    '^[Ee]nthalpy': 'enthalpy',
    '^[d]Density': 'density',
}

varDict_DiffusionEquation = {
    '^[Tt]emperature': 'temperature',
}

varDict_TransportEquation = {
    '^[Vv]oid_*[Ff]raction': 'void_fraction',
    '^[Tt]emperature': 'temperature',
    '^[Ee]nthalpy': 'enthalpy',
    '^[d]Density': 'density',
}

varDict_FiveEqsTwoFluid = {
    '^[Vv]oid_*[Ff]raction': 'void_fraction',
    '^[Pp]ressure': 'pressure',
    '^[Vv]elocity1_*[Xx]': 'gas_velocity_x',
    '^[Vv]elocity1_*[Yy]': 'gas_velocity_y',
    '^[Vv]elocity1_*[zZ]': 'gas_velocity_z',
    '^[Vv]elocity2_*[Xx]': 'liquid_velocity_x',
    '^[Vv]elocity2_*[Yy]': 'liquid_velocity_y',
    '^[Vv]elocity2_*[zZ]': 'liquid_velocity_z',
    '^[Tt]emperature': 'temperature',
}

varDict_IsothermalTwoFluid = {
    '^[Vv]oid_*[Ff]raction': 'void_fraction',
    '^[Pp]ressure': 'pressure',
    '^[Vv]elocity1_*[Xx]': 'gas_velocity_x',
    '^[Vv]elocity1_*[Yy]': 'gas_velocity_y',
    '^[Vv]elocity1_*[zZ]': 'gas_velocity_z',
    '^[Vv]elocity2_*[Xx]': 'liquid_velocity_x',
    '^[Vv]elocity2_*[Yy]': 'liquid_velocity_y',
    '^[Vv]elocity2_*[zZ]': 'liquid_velocity_z',
}
