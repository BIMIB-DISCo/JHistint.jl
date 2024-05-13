### -*- Mode: Julia -*-

### JHistintLogger -- JHistint
### JHistintLogger.jl

### Packages
using Logging
using Dates

### Exported functions
export jh_log_message
export jh_close_logger
export jh_open_logger

"""
    jh_open_logger()

Function to setup and open loggers (console and file io)
"""
function jh_open_logger()
    global io = open(joinpath(DirectoryManager.CONFIG_DIR, "JHistintLog.txt"), "a")
    global logger = SimpleLogger(io)
    global console = ConsoleLogger(stdout)
end

"""
    jh_log_message(level::AbstractString, message::AbstractString)

Function to print log message on console and setup messages to print

# Arguments
- `level::AbstractString` = log level message
        [@debug, @info, @warn, @error]
- `message::AbstractString` = message to print
"""
function jh_log_message(level::AbstractString, message::AbstractString)
    global_logger(console)
    message = string(message, " ", Dates.format(now(), RFC1123Format)) 
    expr = Meta.parse("$(Symbol(level))($(repr(message)))")
    eval(expr)
    global_logger(logger)
    eval(expr)
end

"""
    jh_close_logger()

Function to print  on logfile.txt in CONFIG_DIR and close file io
"""
function jh_close_logger()
    flush(io)
    close(io)
end