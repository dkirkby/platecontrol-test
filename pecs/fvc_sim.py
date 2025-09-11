import pandas
import random
from DOSlib.positioner_index import PositionerIndex

class FVC_proxy_sim:
    
    def __init__(self, max_err = 0.01, error_rate=0.01):
        self.max_err = max_err
        self.error_rate = error_rate #Unused, later to have a chance to "not" match a spot
        self.positions = PositionerIndex
        return

    def send_pm_command(self, command, *args, **kwargs):
        return

    def send_spotmatch_command(self, command, *args, **kwargs):
        return

    def send_fvc_command(self, command, *args, **kwargs):
        return

    def generate_reference(self, exptime=None, match_radius=None, reference_file=None, seqid=None):
        return

    def fvcxy_to_qs(self, fvcxy, seqid=None):
        return

    def qs_to_fvcxy(self, qs, seqid=None):
        return

    def measure(self, expected_positions, seqid=None, exptime=None, match_radius=None, all_fiducials=False, matched_only=True):
        """
        Request a list of center positions
        (Illuminator functionality and fiducial control are outside of the scope of measure())

        seqid = unique sequence id (optional, if none is given it will be generated internally)
        expected_positions = List of dictionaries for each actuator. 
                        Dictionary format:  {"id" : device_id, "q" : float, "s" : float, "flags" : uint}
                                   flags    2 : pinhole center 
                                            4 : fiber center 
                                            8 : fiducial center 
                                           32 : bad fiber or fiducial 
                                    q,s : set to the best guess or zero if the FVC should not match this fiber but just return the position; set the flags to 4+32
                        Notes: additional columns required by the FVC (mag, mess_err) will be added internally and initialized to default values.
                               the indices returns for unmatched fibers will be taken from the list of unmatched fiber indixes but they are scrambled 
                               and the position reported is NOT guaranteed to be to position for the fiber actuator with this index.
        match_radius is the radius in fvc pixels the match_center routine using when attempting to match spots to targets.

        It is also possible to pass the expected positions as a pandas dataframe (with columns id, q, s and flags)

        What it will do:
        This function sends the expected positions to the PlateMaker for conversion to fvc coordinates, loads those into the FVC applciations, calls measure and forwards
        the returned list of center position to the PlateMaker for conversion to focal plane coordinates. The converted positions are turned to the caller.
        The routine manages all data format transformations between the FVC and PlateMaker applications.
        Any exceptions will be passed back to the caller.
        """
        # Basic validation of input data
        assert isinstance(expected_positions, list) or \
             isinstance(expected_positions, pandas.core.frame.DataFrame), \
             'expected_positions argument must be a list of dictionaries or pandas dataframe'
        #assert len(expected_positions) != 0, 'expected_positions argument can not be an empty list or dataframe'
        # create an empty dataframe if list is empty (e.g. for debugging. PM will add fiducials)
        if len(expected_positions) == 0:
            dummy = pandas.DataFrame([{'id':'','q':0.0,'s':0.0,'flags':0}])
            expected_positions = dummy.drop(0)
        if isinstance(expected_positions, list):
            assert isinstance(expected_positions[0], dict),'expected_positions argument must be a list of dictionaries'
            # Convert list of dictionaries to dataframe
            expected_positions = pandas.DataFrame(expected_positions)
        elif isinstance(expected_positions, pandas.core.frame.DataFrame):
            #Assume same format as returned from get_positions in petalApp X1 and X2 must be Q and S
            expected_positions.rename(columns = {'X1':'q'}, inplace=True)
            expected_positions.rename(columns = {'X2':'s'}, inplace=True)
            expected_positions.rename(columns = {'FLAG':'flags'}, inplace=True)
            expected_positions.rename(columns = {'FLAGS':'flags'}, inplace=True) # FLAG -> FLAGS is fixed in some places but not others so do both
            expected_positions.rename(columns = {'DEVICE_ID':'id'}, inplace=True)
            #expected_positions.drop(columns=['PETAL_LOC','DEVICE_LOC'], inplace=True)

        # Add a check that the required keys (id, q, s and flags are given)
        expected_positions = expected_positions[['id', 'q', 's', 'flags']]

        measured_positions = []
        for index, row in expected_positions.iterrows():
            sign = 1
            if random.random() > 0.5:
                sign = -1
            dq = random.random()*self.max_err*sign
            sign = 1
            if random.random() > 0.5:
                sign = -1
            ds = random.random()*self.max_err*sign
            row['flags'] = int(row['flags'])
            row['flags'] |= 1 #set matched
            measured_position = {'id':row['id'], 'q':row['q']+dq,'s':row['s']+ds,'dq':dq,'ds':ds,'flags':row['flags'],'fwhm':random.random(),'mag':random.random(),'peak':random.random()}
            measured_positions.append(measured_position)
        return measured_positions 

    def locate(self, *args, **kwargs):
        """
        Convenience function - same as self.send_fvc_command('locate', *args, **kwargs)
        """
        return self.send_fvc_command('locate', *args, **kwargs)

    def match_centers(self, *args, **kwargs):
        """
        Convenience function - same as self.send_spotmatch_command('match_centers', *args, **kwargs)
        """
        return self.send_spotmatch_command('match_centers', *args, **kwargs)

    def make_targets(self, *args, **kwargs):
        """
        Convenience function - same as self.send_fvc_command('make_targets', *args, **kwargs)
        """
        return self.send_fvc_command('make_targets', *args, **kwargs)

    def calibrate_bias(self, *args, **kwargs):
        """
        Convenience function - same as self.send_fvc_command('calibrate_bias', *args, **kwargs)
        """
        return self.send_fvc_command('calibrate_bias', *args, **kwargs)

    def calibrate_image(self, *args, **kwargs):
        """
        Convenience function - same as self.send_fvc_command('calibrate_image', *args, **kwargs)
        """
        return self.send_fvc_command('calibrate_image', *args, **kwargs)

    def get(self, *args, **kwargs):
        """
        Convenience function - same as self.send_fvc_command('get', *args, **kwargs)
        some commands go to the PlateMaker or spotmatch, others to the FVC.
        """
        if 'instrument' in args:
            return 'sim' #self.send_pm_command('get', *args, **kwargs)
        elif 'match_radius' in args:
            return self.send_spotmatch_command('get', 'match_radius')
        else:
            return self.send_fvc_command('get', *args, **kwargs)

    def set(self, *args, **kwargs):
        """
        Convenience function - same as self.send_fvc_command('det', *args, **kwargs)
        some commands go to the PlateMaker, others to the FVC. This needs to be refined
        """
        if 'instrument' in kwargs:
            return self.send_pm_command('det', *args, **kwargs)
        elif 'match_radius' in kwargs:
            try:
                m = float(kwargs['match_radius'])
                return self.send_spotmatch_command('set', match_radius = m)
            except:
                rstring = 'set: invalid parameter for set match_radius command (%r)' % kwargs['match_radius']
                return 'FAILED: ' + rstring
        else:
            return self.send_fvc_command('set', *args, **kwargs)
