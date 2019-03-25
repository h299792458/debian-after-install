# Debian after install
# Install Intel wifi driver if necessary
sudo apt-get install firmware-iwlwifi
# Config sudo with new file
visudo -f /etc/sudoers.d/myusername
# Add one line
myusername ALL=(ALL) ALL
# If the Xfce desktop enviroment is used, chinese fonts installation is required (synaptic)
# add the following to the .bashrc file if chinese display is desired
export LANG=zh_CN-UTF8
# optional
sudo dpkg-reconfigure locales
# Input method
sudo apt-get install fcitx fcitx-sunpinyin
im-config
# Edit update sources list "/etc/apt/sources.list"
deb http://deb.debian.org/debian/ stretch main contrib non-free
deb-src http://deb.debian.org/debian/ stretch main contrib non-free
deb http://security.debian.org/ stretch/updates main contrib non-free
deb http://deb.debian.org/debian/ stretch-updates main contrib non-free
# Update
sudo apt-get update && sudo apt-get upgrade
# Perform a quick scan for open ports
sudo nmap -sS localhost
# [Option] Enable BBR by adding two lines in /etc/sysctl.conf
net.core.default_qdisc=fq
net.ipv4.tcp_congestion_control=bbr
# Make a firewall
# Install applications
sudo apt-get install texlive
# Install virtualbox vm from offical website
## WARNING: The vboxdrv kernel module is not loaded.
## fix: using synaptic to fix broken packages...
sudo apt-get install linux-headers-$(uname -r)
sudo /sbin/vboxconfig
# Install Matlab
sudo mount -t auto -o loop ~/path-to-matlab/Matlab_2016a/R2016a_glnxa64.iso ~/matlab/
# Matlab desktop entry
[Desktop Entry]
Exec=/usr/local/MATLAB/R2016a/bin/matlab -desktop
Icon=/usr/local/MATLAB/R2016a/toolbox/shared/dastudio/resources/MatlabIcon.png
Terminal=false
StartupNotify=true
