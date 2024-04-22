# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 14:32:19 2024

@author: neuro
"""
# Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch18_ephys\iRNs\GADi39\240206_C+F_2399\Behaviour

import cv2
import os
import pathlib as pl

def avi_to_mp4(input_path, output_path, codec='H264'):
    # Read the .avi file
    cap = cv2.VideoCapture(input_path)

    # Get video properties
    fps = int(cap.get(cv2.CAP_PROP_FPS))
    frame_size = (int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)), int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT)))

    # Define codec and create VideoWriter
    fourcc = cv2.VideoWriter_fourcc(*codec)
    out = cv2.VideoWriter(output_path, fourcc, fps, frame_size)

    while True:
        ret, frame = cap.read()
        if not ret:
            break
        out.write(frame)

    # Release resources
    cap.release()
    out.release()
    cv2.destroyAllWindows()


codec_choice = 'H264'  # Change to 'H265' if desired

videos_path = list(pl.Path('.').glob('**/*.avi'))
videos_path = [v.as_posix() for v in videos_path]

video_out = videos_path.partition('.')

# avi_to_mp4(input_file, output_file, codec=codec_choice)
# print(f"Conversion complete. Output saved as {output_file}")