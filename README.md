# Debian after install
## Install intel wifi driver if necessary
```sh
sudo apt-get install firmware-iwlwifi
```
## Config sudo with new file
```sh
visudo -f /etc/sudoers.d/username
```

Add one line:
```sh
username ALL=(ALL) ALL
```
## For Xfce DE, chinese fonts installation is required (using synaptic suggested).

Add the following to the .bashrc file if chinese display is desired:
```sh
export LANG=zh_CN-UTF8
```
## _optional_
sudo dpkg-reconfigure locales
## Input method
```sh
sudo apt-get install fcitx fcitx-sunpinyin
im-config
```
## Edit "/etc/apt/sources.list"
```sh
deb http://deb.debian.org/debian/ stretch main contrib non-free
deb-src http://deb.debian.org/debian/ stretch main contrib non-free
deb http://security.debian.org/ stretch/updates main contrib non-free
deb http://deb.debian.org/debian/ stretch-updates main contrib non-free
```
## Update
```sh
sudo apt-get update && sudo apt-get upgrade
```
## Perform a quick scan for open ports
```sh
sudo nmap -sS localhost
```
## [Option] Enable BBR by adding two lines in /etc/sysctl.conf
```sh
net.core.default_qdisc=fq
net.ipv4.tcp_congestion_control=bbr
```
## Make a firewall
## Install applications
sudo apt-get install texlive
## Install virtualbox vm from offical website
# WARNING: The vboxdrv kernel module is not loaded.
# fix: using synaptic to fix broken packages...
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
