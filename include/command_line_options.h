#ifndef COMMAND_LINE_OPTIONS_H
#define COMMAND_LINE_OPTIONS_H

#include<string>
#include<vector>
#include<iostream>

// from https://github.com/jamesgregson/utilities
#include<string_utils.h>
#include<command_line_arguments.h>

class command_line_options {
private:
    class option {
    public:
        bool            found;
        bool            required;
        std::string     name;
        arguments_type  type;
        int             count;
        void           *data;
        std::string     description;
        
        option(){
            found = false;
        }
    };
    
    bool                  m_has_optional;
    bool                  m_has_required;
    size_t                m_longest_option;
    size_t                m_longest_type;
    std::vector< option > m_options;
    
    std::string arg_type( const arguments_type type, const int count ){
        switch( type ){
            case ARGUMENT_FLAG:
                return str_format( "" );
                break;
            case ARGUMENT_BOOL:
                return count == 1 ? str_format("<true|false>") : str_format( "<'true'/'false' X%d>", count );
                break;
            case ARGUMENT_INT:
                return count == 1 ? str_format("<Integer>") : str_format( "<Integer X%d>", count );
                break;
            case ARGUMENT_DOUBLE:
                return count == 1 ? str_format("<Real>") : str_format( "<Real X%d>", count );
                break;
            case ARGUMENT_FLOAT:
                return count == 1 ? str_format("<Real>") : str_format( "<Real X%d>", count );
                break;
            case ARGUMENT_STRING:
                return str_format("<String>");
                break;
        }
        return "???";
    }
public:
    command_line_options(){
        m_longest_option = 0;
        m_longest_type = 0;
    }
    void add_optional_parameter( const std::string &name, const arguments_type type, void *data, const int count, const std::string &description ){
        option newopt;
        newopt.name  = name;
        newopt.type  = type;
        newopt.data  = data;
        newopt.count = count;
        newopt.description = description;
        newopt.found = false;
        newopt.required = false;
        m_options.push_back( newopt );
        m_longest_option = std::max( m_longest_option, name.length() );
        m_longest_type = std::max( m_longest_type, arg_type( type, count ).length() );
        m_has_optional = true;
    }
    
    void add_required_parameter( const std::string &name, const arguments_type type, void *data, const int count, const std::string &description ){
        option newopt;
        newopt.name  = name;
        newopt.type  = type;
        newopt.data  = data;
        newopt.count = count;
        newopt.description = description;
        newopt.found = false;
        newopt.required = true;
        m_options.push_back( newopt );
        m_longest_option = std::max( m_longest_option, name.length() );
        m_longest_type = std::max( m_longest_type, arg_type( type, count ).length() );
        m_has_required = true;
    }
    
    bool parse( const int argc, const char **argv ){
        bool ok = true;
        for( int i=0; i<m_options.size(); i++ ){
            if( m_options[i].required ){
                if( command_line_get_argument( argc, argv, m_options[i].name.c_str(), m_options[i].type, m_options[i].data, m_options[i].count ) )
                    m_options[i].found = true;
                else
                    ok = false;
            } else {
                if( command_line_get_argument( argc, argv, m_options[i].name.c_str(), m_options[i].type, m_options[i].data, m_options[i].count ) )
                    m_options[i].found = true;
            }
        }
        if( !ok ){
            std::cout << "==================================================" << std::endl;
            std::cout << "Usage: " << argv[0] << " [parameters]" << std::endl;
            std::cout << "==================================================" << std::endl;
            if( m_has_required || m_has_optional ){
                std::cout << "Parameters:" << std::endl;
                
                if( m_has_required ){
                    for( int i=0; i<m_options.size(); i++ ){
                        if( m_options[i].required ){
                            std::string type = arg_type( m_options[i].type, m_options[i].count );
                            type.resize( m_longest_type, ' ');
                            std::string tmp = m_options[i].name;
                            tmp.resize( m_longest_option, ' ' );
                            std::cout << "  " << tmp << " " << type << " [REQUIRED] " << m_options[i].description << std::endl;
                        }
                    }
                }
                if( m_has_optional ){
                    std::cout << "Optional parameters:" << std::endl;
                    for( int i=0; i<m_options.size(); i++ ){
                        if( !m_options[i].required ){
                            std::string type = arg_type( m_options[i].type, m_options[i].count );
                            type.resize( m_longest_type, ' ');
                            std::string tmp = m_options[i].name;
                            tmp.resize( m_longest_option, ' ' );
                            std::cout << "  " << tmp << " " << type << " [OPTIONAL] " << m_options[i].description << std::endl;
                        }
                    }
                }
            }
            std::cout << "==================================================" << std::endl;
            std::cout << "Errors: " << std::endl;
            for( int i=0; i<m_options.size(); i++ ){
                if( m_options[i].required && !m_options[i].found ){
                    std::cout << "required option " << m_options[i].name << " was not specified!" << std::endl;
                }
            }
            std::cout << std::endl << "command line was: " << std::endl;
            for( int i=0; i<argc; i++ ){
                std::cout << argv[i] << " ";
            }
            std::cout << std::endl << std::endl;
            return false;
        }
        return true;
    }
};

#endif
