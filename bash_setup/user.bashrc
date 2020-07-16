# User .bashrc file.  Should live at ~/.bashrc
# User specific aliases and functions
# NOTE: generally, you don't want to add any export of env variables here.  Use .bash_profile

#Set your prompt.  You mostly have a color supporting xterm these days
if [ "$TERM" == "xterm" -o "$TERM" == "xterm-256color" ]; then
        PS1='\[\033]0;\u@\h: \w\007\]\u@\h[\W]% '
        alias ls='ls --color=tty'
        alias l='ls -aF --color=tty'
        alias ll='ls -lF --color=tty'
else
        PS1='\u@\h[\W]% '
        alias l='ls -aF '
        alias ll='ls -lF'
fi
# Tell bash to use fancy tab completions
if [ -f /etc/bash_completion ]; then
        . /etc/bash_completion
fi
#some useful aliases:
alias h='history | more'
alias j='jobs -l'
alias screen='TERM=screen screen'
alias condorq='condor_q -all -nobatch'
alias cstat='condor_status -constraint "PartitionableSlot =?= TRUE" -format "%32s" Name -format "%6d" TotalSlotCpus -format "%6d" Cpus -format "%5d%%" "((TotalSlotCpus - Cpus) / (TotalSlotCpus * 1.0)) * 10
0" -format "%8d" TotalSlotMemory -format "%8d" Memory -format "%5d%%" "((TotalSlotMemory - Memory) / (TotalSlotMemory * 1.0)) * 100" -format " %4.2f" TotalLoadAvg -format "%3d" DetectedGPUs -format "\n" TR
UE'
ff () { find . -name ${1} -print ; }
