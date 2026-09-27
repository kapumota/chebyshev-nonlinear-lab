function(enable_sanitizers target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} INTERFACE
      -fsanitize=address,undefined -fno-omit-frame-pointer -g)
    target_link_options(${target} INTERFACE -fsanitize=address,undefined)
  endif()
endfunction()
