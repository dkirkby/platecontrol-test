import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import os
import time
import posconstants as pc

class PosAnimator(object):
    """Handles plotting animated visualizations for array of positioners.

    Saving movie files of animations is implemented. Must have ffmpeg to write the
    movie. For Windows machines, download a static binary of ffmpeg, and copy just
    the executable file 'ffmpeg.exe' into the working directory. As of 2016-01-31,
    binaries of ffmpeg are available at: http://ffmpeg.zeranoe.com/builds
    """
    def __init__(self, fignum=0, timestep=0.1):
        self.live_animate = False # whether to plot the animation live
        self.save_movie = True # whether to write out a movie file of animation
        self.ffmpeg_path = os.path.join(pc.petal_directory,'ffmpeg')
        if not (os.path.exists(self.ffmpeg_path + '.exe') or os.path.exists(self.ffmpeg_path)):
            self.ffmpeg_path = 'ffmpeg' # hope the user has it on sys path in this case
        self.codec = 'libx264' # video codec for ffmpeg to use. alt 'libx265'
        self.delete_imgs = False # whether to delete individual image files after generating animation
        self.save_dir = os.path.join(pc.dirs['temp_files'], 'schedule_animations')
        self.frame_dir = '' # generated automatically when saving frames
        self.filename_suffix = '' # optional for user to add before generating animation
        self.add_timestamp_prefix_to_filename = True
        self.framefile_prefix = 'frame'
        self.framefile_extension = '.png'
        self.n_framefile_digits = 5
        self.fignum = fignum
        self.timestep = timestep # frame interval for animations
        self.start_end_still_time = 0.0 # seconds, can add extra still frames at start / end for this duration; gives a visual pause between move sequences
        self.items = {} # keys shall identify the individual items that get drawn, i.e. 'ferrule 321', 'phi arm 42', 'GFA', etc
                        # values are sub-dictionaries, with the entries:
                        #  'time'  : [] # list of time values at which an update or change is to be applied
                        #  'poly'  : [] # list of arrays of polygon points (each is 2xN), defining the item polygon to draw at that time
                        #  'style' : [] # list of dictionaries, defining the plotting style to draw the polygon at that time
        self.global_notes = {0: ''}
        self.labels = {}
        self.label_size = 'x-small'
        self.cropping_on = True # crop the frames to just surround the robots
        self.crop_margin = 14.0 # mm
        self.crop_box = {'xmin':-np.inf, 'xmax':np.inf, 'ymin': -np.inf, 'ymax':np.inf}
        self.pospoly_keys = {'ferrule', 'phi arm', 'central body', 'line at 180', 'Eo', 'Ei', 'Ee'}
        self.fixpoly_keys = {'PTL','GFA'}
        self.styles = pc.plot_styles


    def clear(self):
        """Clear the animator completely of old data.
        """
        self = PosAnimator(self.fignum, self.timestep)
    
    def clear_after(self, time):
        '''Clear the animator of existing data from value time (in seconds)
        onward.'''
        for item in self.items.values():
            while item['time'] and item['time'][-1] >= time:
                for key in ['time', 'poly', 'style']:
                    item[key].pop()
        self.global_notes = {time: note for time, note in self.global_notes.items() if time <= time}

    def is_empty(self):
        """Whether the animator contains any frame data yet."""
        return len(self.items) == 0

    def add_or_change_item(self, item_str, item_idx, time, polygon_points, style_override=''):
        """Add a polygonal item at a particular time to the animation data.
            
            item_str       ... valid options are string keys defined in self.pospoly_keys and self.fixpoly_keys
            item_idx       ... numeric index which gets appended to item_str, in particular to distinguish multiple positioners from each other
            time           ... seconds, time at which this item should be shown
            polygon_points ... as retrieved from PosPoly object
            style_override ... valid options are: '' (no override), 'collision', 'frozen'
        """
        key = str(item_str) + ' ' + str(item_idx)
        if key not in self.items:
            item = {'time':[], 'poly':[], 'style':[], 'is pos poly':item_str in self.pospoly_keys}
        else:
            item = self.items[key]
        if time in item['time']:
            idx = item['time'].index(time)
            replace_existing_timestep = True
        else:
            replace_existing_timestep = False
            if len(item['time']) == 0 or time >= max(item['time']):
                idx = len(item['time'])
            elif time < min(item['time']):
                idx = 0
            else:
                for i in range(len(item['time'])):
                    if item['time'][i] > time:
                        idx = i - 1  # this is where the new timestep data will get inserted
                        break
        if style_override:
            style = self.styles[style_override]
        else:
            style = self.styles[item_str]
        if replace_existing_timestep:
            item['poly'][idx] = polygon_points
            item['style'][idx] = style
        else:
            item['time'].insert(idx, time)
            item['poly'].insert(idx, polygon_points)
            item['style'].insert(idx, style)            
        self.items[key] = item
        
    def add_label(self, text, x, y):
        """Add a text string at position (x,y)."""
        key = len(self.labels)
        self.labels[key] = {'text':text, 'x':x, 'y':y}

    @property    
    def all_times(self):
        """Inspect contents and return the frame times."""
        temp = np.array([])
        for item in self.items.values():
            temp = np.append(temp, item['time'])
        all_times = np.unique(temp)        
        return all_times
    
    def set_note(self, note=None, time=None):
        '''Add a string to the animation plot.
        
            note ... str, will be displayed on plots. None clears existing note.
            time ... seconds, global time. None sets note to display at latest time.
        '''
        note = str(note) if note else ''
        time = float(time) if time else max(self.all_times)
        self.global_notes[time] = note

    def anim_init(self):
        """Sets up the animation, using the data that has been already entered via the
        add_or_change_item method. Returns boolean stating success or not.
        """
        if self.live_animate:
            plt.ion()
        else:
            plt.ioff()
        self.anim_fig = plt.figure(self.fignum, figsize=(20,15))
        self.anim_ax = plt.axes()
        all_times = self.all_times
        if len(all_times) == 0:
            return False
        self.finish_time = max(all_times)
        self.patches = []
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf
        i = 0
        for item in self.items.values():
            patch = self.get_patch(item,0)
            self.patches.append(self.anim_ax.add_patch(patch))
            if (self.cropping_on and item['is pos poly']) or not self.cropping_on:
                margin = self.crop_margin if self.cropping_on else 0.0
                xmin = min(xmin,min(item['poly'][0][0]) - margin)
                xmax = max(xmax,max(item['poly'][0][0]) + margin)
                ymin = min(ymin,min(item['poly'][0][1]) - margin)
                ymax = max(ymax,max(item['poly'][0][1]) + margin)
            item['patch_idx'] = i
            item['last_frame'] = all_times.tolist().index(item['time'][-1])
            i += 1
        for label in self.labels.values():
            plt.text(s=label['text'], x=label['x'], y=label['y'], family='monospace', horizontalalignment='center', size=self.label_size)
        plt.axis('square')
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        return True

    def anim_frame_update(self, frame):
        """Sequentially update the animation to the next frame.
        """
        for item in self.items.values():
            i = item['patch_idx']
            if frame <= item['last_frame'] and len(item['poly']) > frame:
                self.set_patch(self.patches[i], item, frame)

    def animate(self):
        '''Returns path of output file.'''
        successful = self.anim_init()
        if not successful:
            print('Animator not initialized. Usually due to no frames available to animate.')
            return
        frame_number = 1
        stdout_message_period = 50 # number of frames per update message
        image_paths = {}
        if self.add_timestamp_prefix_to_filename:
            timestamp = pc.filename_timestamp_str() + '_'
        else:
            timestamp = ''
        suffix = '_' + str(self.filename_suffix) if self.filename_suffix else ''
        self.frame_dir = os.path.join(self.save_dir, timestamp + 'frames' + suffix)
        if self.live_animate:
            plt.show()
        if self.save_movie:
            fps = 1/self.timestep
            if not(os.path.exists(self.save_dir)):
                os.mkdir(self.save_dir)
            if not(os.path.exists(self.frame_dir)):
                os.mkdir(self.frame_dir)
            for i in range(round(self.start_end_still_time/self.timestep)):
                frame_number,path = self.grab_frame(frame_number)
                image_paths[frame_number] = path
        all_times = self.all_times
        frame_times = np.arange(min(all_times), max(all_times)+self.timestep/2, self.timestep) # extra half-timestep on max ensures inclusion of max val in range
        note_times = sorted(self.global_notes.keys())
        assert max(note_times) <= max(frame_times)
        assert min(note_times) >= min(frame_times)
        note_idx = 0
        for frame in range(len(frame_times)):
            frame_start_time = time.perf_counter()
            self.anim_frame_update(frame)
            if note_times:
                this_time = frame_times[frame]
                until_now = [i for i in range(len(note_times)) if note_times[i] <= this_time]
                note_idx = until_now[-1]
                note = self.global_notes[note_times[note_idx]]
            else:
                note = ''
            title = f'time: {frame*self.timestep:5.2f} / {self.finish_time:5.2f} sec'
            if note:
                title += f'\n{note}'
            plt.title(title)
            if self.save_movie:
                frame_number,path = self.grab_frame(frame_number)
                image_paths[frame_number] = path
                if frame % stdout_message_period == 0:
                    print(' ... animation frame ' + str(frame) + ' of approx ' + str(len(frame_times)) + ' saved')
            if self.live_animate:
                plt.pause(0.001) # just to force a refresh of the plot window (hacky, yep. that's matplotlib for ya.)
                time.sleep(max(0,self.timestep-(time.perf_counter()-frame_start_time))) # likely ineffectual (since matplotlib so slow anyway) attempt to roughly take out frame update / write time
        if self.save_movie:
            for i in range(round(self.start_end_still_time/self.timestep)):
                frame_number,path = self.grab_frame(frame_number)
                image_paths[frame_number] = path
        plt.close()
        if self.save_movie:
            input_file = os.path.join(self.frame_dir, self.framefile_prefix + '%' + str(self.n_framefile_digits) + 'd' + self.framefile_extension)
            output_file = os.path.join(self.save_dir, timestamp + 'schedule_anim' + suffix + '.mp4')
            ffmpeg_cmd = self.ffmpeg_path + ' -y -r ' + str(fps) + ' -i ' + input_file + ' -vcodec ' + self.codec + ' ' + output_file
            err = os.system(ffmpeg_cmd)
            if err:
                output_file = 'FAILED - check ffmpeg installation'
        if self.delete_imgs:
            for path in image_paths.values():
                os.remove(path)
        return output_file

    def grab_frame(self, frame_number):
        """Saves current figure to an image file.
        Returns a tuple containing:
            next frame number
            path saved to
        """
        path = os.path.join(self.frame_dir, self.framefile_prefix + str(frame_number).zfill(self.n_framefile_digits) + self.framefile_extension)
        plt.savefig(path, bbox_inches='tight')
        return frame_number + 1, path  
    
    @staticmethod
    def get_patch(item,index):
        return plt.Polygon(pc.transpose(item['poly'][index]),
                           linestyle=item['style'][index]['linestyle'],
                           linewidth=item['style'][index]['linewidth'],
                           edgecolor=item['style'][index]['edgecolor'],
                           facecolor=item['style'][index]['facecolor'])

    @staticmethod
    def set_patch(patch,item,index):
        patch.set_xy(pc.transpose(item['poly'][index]))
        patch.set_linestyle(item['style'][index]['linestyle'])
        patch.set_linewidth(item['style'][index]['linewidth'])
        patch.set_edgecolor(item['style'][index]['edgecolor'])
        patch.set_facecolor(item['style'][index]['facecolor'])
