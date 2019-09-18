# Debian after install


Several things to do after installing debian


## Edit "/etc/apt/sources.list"
Add sources
```sh
deb http://mirrors.ustc.edu.cn/debian buster main contrib non-free
deb http://security.debian.org/debian-security/ buster/updates main contrib non-free
```
Update and upgrade
```sh
sudo apt-get update && sudo apt-get upgrade
```


## Config sudo
Create a new config file
```sh
visudo -f /etc/sudoers.d/username
```
Add one line
```sh
username ALL=(ALL) ALL
```


## Check the OS release information
```sh
cat /etc/os-release
```


## Check the CPU info
```sh
sudo less /proc/cpuinfo
```


## Update CPU microcode
```sh
sudo dmesg | grep microcode
sudo apt-get install intel-microcode
sudo reboot
```


## Install Intel wifi driver
```sh
apt-get install firmware-iwlwifi
modprobe -r iwlwifi && modprobe iwlwifi
```


## 中文字体
Choose and install free Chinese fonts
```sh
apt-get instal fonts-wqy-zenhei
```
1、进入windows系统，到C:\Windows\Fonts目录，复制拷贝自己需要的字体（也可以全部拷贝，包含windows支持的所有中文字体）;
2、在debian系统下，将需要的字体文件拷贝到libreoffice的配置文件路径～/.config/libreoffice/4/user/fonts下。如果没有fonts目录新建即可。

Add the following to the .bashrc file if Chinese display is desired
```sh
export LANG=zh_CN-UTF8
```
_optional_
```sh
sudo dpkg-reconfigure locales
```


## 中文输入法
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
sudo apt-get install texlive-pictures
sudo apt-get install texlive-science
sudo apt-get install texlive-latex-extra
```


## Install virtualbox
Download from the offical website suggested

How to fix
```sh
WARNING: The vboxdrv kernel module is not loaded.
```
run the following
```sh
sudo apt-get install linux-headers-$(uname -r)
sudo /sbin/vboxconfig
```


## Mount NTFS
```sh
sudo apt-get install ntfs-3g
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

