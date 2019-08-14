# Debian after install

Several things to do after installing debian

## Check the OS release information
```sh
cat /etc/os-release
```

## Edit "/etc/apt/sources.list"
Add sources
```sh
deb http://mirrors.ustc.edu.cn/debian buster main contrib non-free
deb http://security.debian.org/debian-security/ buster/updates main contrib non-free
```
Update and upgrade
```sh
apt-get update && apt-get upgrade
```

## Intel wifi driver
```sh
apt-get install firmware-iwlwifi
modprobe -r iwlwifi
modprobe iwlwifi
```

## Install sudo
Create a new config file
```sh
visudo -f /etc/sudoers.d/username
```
Add one line
```sh
username ALL=(ALL) ALL
```

## Install language support
Choose and install Chinese fonts
```sh
apt-get instal fonts-wqy-zenhei
```

Add the following to the .bashrc file if Chinese display is desired
```sh
export LANG=zh_CN-UTF8
```

_optional_
```sh
sudo dpkg-reconfigure locales
```

## Input method
```sh
sudo apt-get install fcitx fcitx-sunpinyin
```

_optional_
```sh
im-config
```

## Network
Perform a quick scan for open ports:
```sh
sudo nmap -sS localhost
```

_optional_
Enable BBR by adding two lines in /etc/sysctl.conf:
```sh
net.core.default_qdisc=fq
net.ipv4.tcp_congestion_control=bbr
```

## Make a firewall
TODO

## Install Latex
```sh
sudo apt-get install texlive
```

## Install virtualbox
Download from the offical website suggested

How to fix
```sh
WARNING: The vboxdrv kernel module is not loaded.
```
1.  fix bad package independences via e.g. Synaptic.
2.  run the following:
```sh
sudo apt-get install linux-headers-$(uname -r)
sudo /sbin/vboxconfig
```

## Install Matlab
Mount the iso file:
```sh
sudo mount -t auto -o loop ~/path-to-matlab/Matlab_2016a/R2016a_glnxa64.iso ~/matlab/
```

Create desktop entry:
```sh
[Desktop Entry]
Exec=/usr/local/MATLAB/R2016a/bin/matlab -desktop
Icon=/usr/local/MATLAB/R2016a/toolbox/shared/dastudio/resources/MatlabIcon.png
Terminal=false
StartupNotify=true
```

## Mount Android phone (MTP)
mount:
```sh
jmtpfs Android/
```
unmount:
```sh
fusermount -u Android/
```
