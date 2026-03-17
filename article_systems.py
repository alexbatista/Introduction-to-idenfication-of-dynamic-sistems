import control as ct
import pandas as pd

s = ct.tf('s') #ct.TransferFunction.s

# Process model, g0(s) E1, E4, E7, E10, E13 given in the article

E1 = 1/((s+1)*(0.2*s+1))

E4 = 1/(s+1)**4

E7 = (-2*s + 1)/(s+1)**3

delay = ct.tf(*ct.pade(1, 5))  # 5th-order Padé approximation of e^{-s}

E10 = delay / ((20*s+1)*(2*s+1))

E13 = ((2*s+1) / ((10*s+1)*(0.5*s+1))) * delay

functions = dict(e1=E1,e4=E4,e7=E7,e10=E10,e13=E13)

for key,value in functions.items():
    print(f"Transfer function G(s) of {key.upper()}: ")
    print(value)



for key,value in functions.items():
    time, response = ct.step_response(value)
    df = pd.DataFrame({'Time': time, 'Output': response})
    df.to_csv(f'dataset/step_response_{key.upper()}.csv', index=False, header=False)
    print(f"Step response data saved to dataset folder as step_response_{key.upper()}")

