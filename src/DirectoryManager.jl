### -*- Mode: Julia -*-

### DirectoryManager -- JHistint
### DirectoryManager.jl

### Includes the logic to save and update the current user workspace.

module DirectoryManager

### Packages

using JSON


### Exported functions

export set_environment


### Exported constants

export CONFIG_DIR, DATA_DIR

### Identify the standard directories for data and resources based on
### the operating system.

if Sys.islinux()
    
    ## According to XDG Standards: $XDR_CONFIG_HOME
    
    const CONFIG_DIR = joinpath(homedir(), ".config", "JHistint")
    
    ## According to XDG Standards: $XDR_DATA_HOME is $HOME/.local/share

    ## Initialize the standard Workspace directory if none has been set
    const DATA_DIR = joinpath(homedir(), ".local", "share", "JHistint")

elseif Sys.iswindows()

    ## Config data will be saved in %APPDATA% into Local folder
    
    const CONFIG_DIR = joinpath(homedir(), "AppData", "Local", "JHistint")
    const DATA_DIR = CONFIG_DIR
    
elseif Sys.isapple()
    const CONFIG_DIR = joinpath(homedir(), ".config", "JHistint")
    const DATA_DIR = CONFIG_DIR

else
    ## For all the other operating systems the configuration files
    ## will be saved in JHisting.jl folder

    const CONFIG_DIR = @__DIR__
    const DATA_DIR = joinpath(homedir(), "JHistint")
end


### Initialize the standard directories on first launch of the app on
### the system.

if !isdir(CONFIG_DIR) 
    mkdir(CONFIG_DIR)
end

if !isdir(DATA_DIR) 
    mkdir(DATA_DIR)
end


function set_environment()
    if !isdir(joinpath(CONFIG_DIR, "collections"))
        mkdir(joinpath(CONFIG_DIR, "collections"))
    end
    if !isdir(joinpath(CONFIG_DIR, "cases"))
        mkdir(joinpath(CONFIG_DIR, "cases"))
    end
end

end # module JHistint.Workspace

### DirectoryManager.jl ends here.
