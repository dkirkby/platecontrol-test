This folder contains software to communicate with Mightex LED Controllers via USB on Linux systems.
The supplied executable, mightex_cmd, was compiled for Lubuntu.
No attempt has yet been made to compile on other versions of Linux, but I believe the code is fairly
generic, so there's a decent chance it will compile and work.

mightex_cmd communicates with the LED controller via the /dev/hidraw* devices.
These devices are generally not readable and writeable by normal users.
To overcome that, place the 19-mightex_udev.rules file into /etc/udev/rules.d
as root (sudo). This will make the hidraw Mightex devices readable and writeable by everyone.

The file mightex_cmd.help.txt contains usage instructions and some examples for the mightex_cmd program.

From Python3, typical useage will be via the class defined for Mightex LED Controllers in mightex.py.
NOTE: mightex.py assumes that mightex_cmd is in your PATH!

Here are the methods currently defined for the class Mightex_LED_Controller.
Note: failures generally raise IOError.

    # Open a specific channel (default 1) on the Mightex LED Controller with the specified
    # serial number, or if no serial number is specified, on the last LED controller in the list
    # (i.e., essentially chosen at random since the order on the list is determined by the
    # order of connection/discovery)
    def __init__(self,channel=1,serialno=""):

    # returns the maximum number of channels for the current Mightex LED Controller
    def getMaxChannels(self):

    # returns a list of the serial numbers for currently attached Mightex LED Controllers
    def getSerialNos(self):

    # returns the present mode of the current Channel
    def getMode(self):

    # set the mode. returns 1 on success, 0 on failure. Mode=0 => Off, Mode=1 => Normal
    def setMode(self,mode=0):

    # returns a list of the Maximum and Current milliAmp levels for the current Channel        
    def getLevels(self):

    # Set the value of the LED current in milliAmps
    # Optionally, specify the Mode and the Maximum LED current
    # If the Mode is not specified, none will be set (NOTE: no changes to the
    # LED intensity are made until you set Mode 1, even if already in Mode 1)
    # If the Maximum is not specified, it is set to 10 more than what you specified for the LED Current
    # (up to a maximum of 1000)
    def setLevel(self,led_milliamps=0,mode=-1,max_led_milliamps=-1):

    # Load the Factory Defaults for all settings
    # NOTE: the new settings do not take effect until you set the Mode
    # and they are NOT saved in the NVRAM
    def setFactoryDefaults(self):

    # Reset the Mightex LED Controller
    # CAUTION: Executing this command may require a Linux restart to recover!
    def setReset(self):

Typical usage of the Python3 class might be something like this:
import mightex

m=mightex.Mightex_LED_Controller(1,'02-170175-002') # open channel 1 of SN="02-170175-002"
m.setLevel(70,1)    # set channel 1 to 70 milliAmps

