import streamlit as st
import inspect

import FundamentalConcepts
import SolidMaterials
import SurfaceRoughness
import NonConformalSolidContact
import LiquidProperties
import LubricantCharacterization
import LubricationConformalContacts
import LubricationNonConformalContacts
import DryAndMixedLubrication
import Wear
import MeasuringFrictionAndWear

functions_dict = {

    "Fundamental Concepts": {
        "Friction": FundamentalConcepts.friction,
        "Lambda Ratio": FundamentalConcepts.lambda_ratio,
        "Hersey Number": FundamentalConcepts.hersey_number
    },

    "Solid Materials": {
        "Shear Modulus": SolidMaterials.shear_modulus,
        "Shear Strength": SolidMaterials.shear_strenth,
        "Hardness": SolidMaterials.hardness
        },

    "Surface Roughness": {
        "Gaussian Distribution": SurfaceRoughness.gaussian_dist,
        "Zi' Mean": SurfaceRoughness.zi_prime_mean,
        "Zi": SurfaceRoughness.zi,
        "Ra": SurfaceRoughness.Ra,
        "Rq": SurfaceRoughness.Rq,
        "Sigma_s": SurfaceRoughness.sig_s,
        "Ra Composite": SurfaceRoughness.Ra_composite,
        "Rq Composite": SurfaceRoughness.Rq_composite,
        "Rsk - Skewness": SurfaceRoughness.Rsk,
        "Rku - Kurtosis": SurfaceRoughness.Rku,
        "Autocorrelations": SurfaceRoughness.autocorrelation,
        "Exponential Autocorrelation": SurfaceRoughness.exp_autocorrelation,
        "Orientation": SurfaceRoughness.orientation
    },
    "Non-Conformal Solid Contact": {
        "Area of a journal bearing": NonConformalSolidContact.area_journal_bearing,
        "Applied Pressure": NonConformalSolidContact.applied_pressure,
        "Effective radius": NonConformalSolidContact.effective_radius,
        "Effective elastic modulus": NonConformalSolidContact.effective_elastic_modulus,
        "a (Point Contact)": NonConformalSolidContact.point.a,
        "Max pressure (Point Contact)": NonConformalSolidContact.point.pmax,
        "Average pressure (Point Contact)": NonConformalSolidContact.point.pavg,
        "Max subsurface shear stress": NonConformalSolidContact.point.tmax,
        "Zmax position (Point Contact)": NonConformalSolidContact.point.zmax,
        "b (Line Contact)": NonConformalSolidContact.line.b,
        "Max pressure (Point Contact)": NonConformalSolidContact.line.pmax,
        "Average pressure (Point Contact)": NonConformalSolidContact.line.pavg,
        "Max subsurface shear stress": NonConformalSolidContact.line.tmax,
        "Zmax position (Point Contact)": NonConformalSolidContact.line.zmax,
        "Plastic area": NonConformalSolidContact.Aplastic,
        "Adhesive force": NonConformalSolidContact.Wadh,
        "Real contact area": NonConformalSolidContact.Ar,
        "Smooth surface approx.": NonConformalSolidContact.x_smooth_surface_approx,
        "Plasticity index": NonConformalSolidContact.plasticity_index
    },
    "Liquid Properties": {
        "Absolute Viscosity": LiquidProperties.absolute_visc,
        "Kinematic Viscosity": LiquidProperties.kinematic_viscosity_v,
        "Williams-Landel-Ferry Equation": LiquidProperties.WilliamsLandelFerry,
        "Vogel Equation": LiquidProperties.Vogel,
        "Barus Equation": LiquidProperties.Barus,
        "Roelands Equation": LiquidProperties.Roelands,
        "Roelands NPT": LiquidProperties.Roelands_npt,
        "PSSI": LiquidProperties.PSSI,
        "Carreau Equation": LiquidProperties.Carreau,
        "Specific Gravity": LiquidProperties.SpecificGravity,
        "Density with temperature": LiquidProperties.DensityWithTemp,
        "Dowson-Higginson Equation": LiquidProperties.DowsonHigginson,
        "Specific Heat": LiquidProperties.SpecificHeat,
        "Thermal Conductivity": LiquidProperties.ThermalConductivity
    },
    "Lubricant Charcaterization": {
        "Viscosity Index": LubricantCharacterization.ViscosityIndex
    },
    "Lubrication of Conformal Contacts": {
        "Film Thickness on Inclined Plane": LubricationConformalContacts.FilmThicknessInclinedPlane,
        "Reynolds - Integrated": LubricationConformalContacts.ReynoldsIntegrated,
        "Max Pressure": LubricationConformalContacts.pmax,
        "Load-Carrying Capacity - Inclined Plane": LubricationConformalContacts.WPrimeInclinedPlane,
        "Viscous Friction - Inclined Plane": LubricationConformalContacts.FvPrimeInclinedPlane,
        "Viscous Friction Ratio": LubricationConformalContacts.fv,
        "Sommerfield Number": LubricationConformalContacts.SommerfieldNumber,
        "pa for Sommerfield": LubricationConformalContacts.pa
    },
    "Lubrication of Non-Conformal Contacts": {
        "SRR": LubricationNonConformalContacts.SRR,
        "Ratio of Radii": LubricationNonConformalContacts.RatioOfRadii,
        "gV": LubricationNonConformalContacts.gV,
        "gE": LubricationNonConformalContacts.gE,
        "h0'": LubricationNonConformalContacts.h0Prime,
        "h0' PE": LubricationNonConformalContacts.h0PrimePE,
        "h0' IR": LubricationNonConformalContacts.h0PrimeIR,
        "h0' PR": LubricationNonConformalContacts.h0PrimePR,
        "h0' IE": LubricationNonConformalContacts.h0PrimeIE
    },
    "Dry and Mixed Friction": {
        "fc": DryAndMixedLubrication.fc,
        "fc Sum": DryAndMixedLubrication.fc_sum,
        "fa Normal": DryAndMixedLubrication.fa_normal,
        "fa Shear Strength": DryAndMixedLubrication.fa_SS,
        "fa Plastic": DryAndMixedLubrication.fa_plastic,
        "fa Elastic": DryAndMixedLubrication.fa_elastic,
        "Deformation Force - Fd": DryAndMixedLubrication.Fd,
        "Deformation Force Ratio fd": DryAndMixedLubrication.fd,
        "fd Ceramic": DryAndMixedLubrication.fd_ceramic,
        "fm": DryAndMixedLubrication.fm,
        "Load Support ξ": DryAndMixedLubrication.xi,
        "Peclet Number": DryAndMixedLubrication.Pe,
        "Average Flash Temp": DryAndMixedLubrication.AverageFlashTemp,
        "Rolling Friction Coefficient fr": DryAndMixedLubrication.fr
    },
    "Wear": {
        "Archard Equation": Wear.Archard,
        "Wi": Wear.Wi,
        "Li": Wear.Li,
        "Vi": Wear.Vi,
        "V": Wear.V,
        "Wi Ductile": Wear.Wi_ductile,
        "Vi Ductile": Wear.Vi_ductile,
        "V Ductile": Wear.V_ductile,
        "V Brittle": Wear.V_brittle
    },
    "Measuring Friction and Wear": {
        "Static Force Coefficient fstatic": MeasuringFrictionAndWear.fstatic
    }
}

def parse_docstring(docstring):
    """Parse the docstring to extract parameter names and descriptions."""
    lines = docstring.strip().split('\n')
    param_info = {}
    instructions = []

    # Check if the first non-empty line is the dashed line
    first_line_is_dashes = True
    for line in lines:
        stripped_line = line.strip()
        if stripped_line == "----------------------------------------------------------------------------------------------------":
            if first_line_is_dashes:
                continue  # Skip this line only if it's the first non-empty line
        else:
            first_line_is_dashes = False  # After seeing a valid line, set this to False

        if stripped_line == "":
            continue
        if ':' in stripped_line:
            param_name, description = stripped_line.split(':', 1)
            param_info[param_name.strip()] = description.strip()
        else:
            instructions.append(stripped_line)

    return instructions, param_info

def generate_inputs(func, tab_name, func_name):
    """Generate input fields based on function parameters and their descriptions."""
    signature = inspect.signature(func)
    instructions, param_info = parse_docstring(func.__doc__)
    inputs = {}

    # Display the instructions in the GUI
    if instructions:
        # st.markdown("### Instructions")
        st.write("\n".join(instructions))

    for param in signature.parameters.values():
        param_name = param.name
        param_type = param.annotation
        description = param_info.get(param_name, "No description available.")

        # Generate input based on parameter type
        if param_type == float:
            default_value = st.session_state.get(f"{tab_name}_{func_name}_{param_name}", 0.0)
            inputs[param_name] = st.number_input(f"{param_name} ({description})", value=default_value, key=f"{tab_name}_{func_name}_{param_name}")
        elif param_type == list:
            default_value = st.session_state.get(f"{tab_name}_{func_name}_{param_name}", "")
            inputs[param_name] = st.text_input(f"{param_name} ({description}) - Enter values separated by commas", value=default_value, key=f"{tab_name}_{func_name}_{param_name}")
        else:
            st.warning(f"Unsupported parameter type: {param_type} for {param_name}")

    return inputs

def clear_inputs(tab_name, func_name):
    """Clear inputs for a specific function."""
    for param in functions_dict[tab_name][func_name].__code__.co_varnames:
        key = f"{tab_name}_{func_name}_{param}"
        if key in st.session_state:
            del st.session_state[key]  # Remove the key from session state

################################################################################

st.title("TriboPy")

# Create tabs for each file
tab_names = list(functions_dict.keys())
# tabs = st.tabs(tab_names)

# Create columns for left and right navigation buttons
col1, col2, col3 = st.columns([1, 8, 1])  # Adjust column sizes as necessary

with col1:
    if st.button("←"):
        # Logic for scrolling left
        if 'tab_index' not in st.session_state:
            st.session_state.tab_index = 0
        else:
            st.session_state.tab_index = max(0, st.session_state.tab_index - 1)

with col2:
    # Show the current tab based on the index
    current_tab = tab_names[st.session_state.get('tab_index', 0)]
    selected_tab = st.selectbox("Select a Topic", options=tab_names, index=st.session_state.get('tab_index', 0), key='tab_selector', on_change=lambda: st.session_state.update({'tab_index': tab_names.index(st.session_state.tab_selector)}))

with col3:
    if st.button("→"):
        # Logic for scrolling right
        if 'tab_index' not in st.session_state:
            st.session_state.tab_index = 0
        else:
            st.session_state.tab_index = min(len(tab_names) - 1, st.session_state.tab_index + 1)

# Iterate through each tab and create UI for each function
for tab_name in tab_names:
    if tab_name == current_tab:
        with st.container():
            st.header(tab_name)  # Use the tab name directly
            
            # Iterate through functions in the current tab
            for func_name, func in functions_dict[tab_name].items():
                st.subheader(func_name)

                # Generate inputs for the function with unique keys
                inputs = generate_inputs(func, tab_name, func_name)

                # Generate a unique key for the buttons
                button_key = f"{tab_name}_{func_name}_calculate"
                clear_key = f"{tab_name}_{func_name}_clear"

                # Calculate button with a unique key
                if st.button("Calculate", key=button_key):
                    try:
                        # Convert list input from text to a list of floats
                        if 'zi' in inputs:
                            inputs['zi'] = [float(z) for z in inputs['zi'].split(',')]
                        
                        # Call the function with unpacked inputs
                        result = func(**inputs)
                        st.success(f"Result: {result}")
                        # Store the result in session state
                        st.session_state[f"{tab_name}_{func_name}_result"] = result
                    except Exception as e:
                        st.error(f"Error: {e}. Please check your input.")

                # Clear button with a unique key
                if st.button("Clear", key=clear_key):
                    clear_inputs(tab_name, func_name)
                    st.session_state[f"{tab_name}_{func_name}_result"] = None  # Clear result display

                # Display the result if it exists
                if f"{tab_name}_{func_name}_result" in st.session_state and st.session_state[f"{tab_name}_{func_name}_result"] is not None:
                    st.info(f"Last Result: {st.session_state[f'{tab_name}_{func_name}_result']}")