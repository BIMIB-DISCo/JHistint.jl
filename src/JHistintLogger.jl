### -*- Mode: Julia -*-

### JHistintLogger -- JHistint
### JHistintLogger.jl

### Exported functions
export info_log
export close_logger
export open_logger

"""
    open_logger()

Function to setup and open loggers (console and file io)
"""
function open_logger()
    global io = open(joinpath(DirectoryManager.CONFIG_DIR, "JHistintLog.txt"), "w+")
    global logger = SimpleLogger(io)
    global console = ConsoleLogger(stdout)
end

"""
    log_message(level::AbstractString, message::AbstractString)

Function to print log message on console and setup messages to print

# Arguments
- `level::AbstractString` = log level message
        [@debug, @info, @warn, @error]
- `message::AbstractString` = message to print
"""
function log_message(level::AbstractString, message::AbstractString)
    global_logger(console)
    message = string(message, " ", Dates.format(now(), RFC1123Format)) 
    expr = Meta.parse("$(Symbol(level))($(repr(message)))")
    eval(expr)
    global_logger(logger)
    eval(expr)
end

"""
function info_log(message::AbstractString)
    global_logger(console)
    message = message, Dates.format(now(), RFC1123Format)
    @info message
    global_logger(logger)
    @info message
end

function warn_log(message::AbstractString)
    global_logger(console)
    message = message, Dates.format(now(), RFC1123Format)
    @warn message
    global_logger(logger)
    @warn message
end

function debug_log(message::AbstractString)
    global_logger(console)
    message = message, Dates.format(now(), RFC1123Format)
    @debug message
    global_logger(logger)
    @debug message
end

function error_log(message::AbstractString)
    global_logger(console)
    message = message, Dates.format(now(), RFC1123Format)
    @error message
    global_logger(logger)
    @error message
end
"""

"""
    close_logger()

Function to print  on logfile.txt in CONFIG_DIR and close file io
"""
function close_logger()
    flush(io)
    close(io)
end