import cv2
import pandas as pd
import os
from tqdm import tqdm
import matplotlib.pyplot as plt

# Open the video file
video_file = r'C:\Users\dmagyar\Desktop\BAfc_videos_pupil-DM-2024-10-25\videos\MD278.avi'  # Replace with your video file
cap = cv2.VideoCapture(video_file)

# Check if the video opened successfully
if not cap.isOpened():
    print("Error: Could not open video.")
    exit()

# Read the first frame to let the user define the ROI
ret, first_frame = cap.read()
if not ret:
    print("Error: Could not read the first frame.")
    exit()

# Allow user to select ROI
x, y, w, h = cv2.selectROI("Select ROI", first_frame, fromCenter=False, showCrosshair=True)

# Close the ROI selection window
cv2.destroyWindow("Select ROI")

# Initialize a list to store intensity values
intensity_values = []

# Get the total number of frames in the video
total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

# Ask the user for the percentage of frames to process
while True:
    try:
        percentage = float(input("Enter the percentage of frames to process (0-100): "))
        if 0 <= percentage <= 100:
            break
        else:
            print("Please enter a valid percentage between 0 and 100.")
    except ValueError:
        print("Invalid input. Please enter a numeric value.")

frames_to_process = int((percentage / 100) * total_frames)

# Create a progress bar using tqdm
print(f"Processing {frames_to_process} frames out of {total_frames} total frames...")
for frame_idx in tqdm(range(frames_to_process), desc="Processing frames", unit="frame"):
    # Read the next frame
    ret, frame = cap.read()
    if not ret:
        break  # Exit loop if no more frames

    # Crop the frame to the ROI
    roi = frame[y:y+h, x:x+w]

    # Calculate total intensity
    total_intensity = roi.sum()  # Sum pixel values
    intensity_values.append(total_intensity)

cap.release()

# Save the intensity values to a CSV file in the same folder as the video
base_filename = os.path.splitext(os.path.basename(video_file))[0]
csv_filename = f"{base_filename}_intensity_values.csv"
output_dir = os.path.dirname(video_file)
csv_path = os.path.join(output_dir, csv_filename)

# Convert the intensity values to a DataFrame and save as CSV
df = pd.DataFrame(intensity_values, columns=['Total Intensity'])
df.to_csv(csv_path, index=False)

print(f'Total intensity values saved to: {csv_path}')

# Calculate mean and standard deviation of the intensity values
mean_intensity = df['Total Intensity'].mean()
std_intensity = df['Total Intensity'].std()
plus_5_sd = mean_intensity + 5 * std_intensity

# Plot the intensity values over time
plt.figure(figsize=(10, 6))
plt.plot(range(len(intensity_values)), intensity_values, label='Total Intensity', color='blue')
plt.axhline(y=plus_5_sd, color='red', linestyle='--', label='+5 SD')

plt.title('Intensity Values Over Time')
plt.xlabel('Frame Index')
plt.ylabel('Total Intensity')
plt.grid()
plt.legend()
plt.show()

# Find frame indices where the intensity crosses the +5 SD threshold
crossing_indices = []
above_threshold = False  # Track if the previous frame was above threshold

for i, intensity in enumerate(intensity_values):
    if intensity > plus_5_sd:
        if not above_threshold:  # Only record if it's the first in a consecutive sequence
            crossing_indices.append(i)
        above_threshold = True
    else:
        above_threshold = False  # Reset if current frame is below threshold

# Create a DataFrame for the crossing indices
crossing_df = pd.DataFrame(crossing_indices, columns=['Frame Index'])

# Save the crossing indices to a new CSV file
crossing_csv_filename = f"{base_filename}_crossing_indices.csv"
crossing_csv_path = os.path.join(output_dir, crossing_csv_filename)
crossing_df.to_csv(crossing_csv_path, index=False)

print(f'Frame indices where intensity crosses +5 SD saved to: {crossing_csv_path}')