# Function to apply Automatic Gain Control (AGC) to the data
def apply_agc(data, window_length=100):
    # Compute the envelope of the data
    envelope = np.abs(data)
    
    # Apply a moving average with the specified window length
    agc_gain = np.convolve(envelope, np.ones(window_length)/window_length, mode='same')
    
    # Avoid division by zero
    agc_gain[agc_gain == 0] = 1
    
    # Apply the gain to the original data
    agc_data = data / agc_gain
    
    return agc_data

# Apply AGC to the transposed data
agc_data = np.apply_along_axis(apply_agc, 1, transposed_data)

# Visualize the AGC-applied data as a radargram using imshow
plt.figure(figsize=(12, 5))
plt.imshow(agc_data, cmap='gray', aspect='auto')
plt.colorbar(label='Amplitude (AGC applied)')
plt.title('AGC-applied Radargram Visualization of 3D_BSCAN_GPR.txt')
plt.xlabel('Trace Number')
plt.ylabel('Depth / Time')
plt.show()


