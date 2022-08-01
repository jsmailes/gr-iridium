/*
 * Copyright 2022 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(header_extract.h)                                          */
/* BINDTOOL_HEADER_FILE_HASH(6b48b04f6cfb0cf139bdd4658e750da8)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <iridium/header_extract.h>
// pydoc.h is automatically generated in the build directory
#include <header_extract_pydoc.h>

void bind_header_extract(py::module& m)
{

    using header_extract    = ::gr::iridium::header_extract;


    py::class_<header_extract, gr::block, gr::basic_block,
        std::shared_ptr<header_extract>>(m, "header_extract", D(header_extract))

        .def(py::init(&header_extract::make),
           py::arg("correlation_sample_rate"),
           py::arg("output_samples"),
           py::arg("search_depth"),
           py::arg("hard_max_queue_len"),
           py::arg("input_taps"),
           py::arg("start_finder_taps"),
           py::arg("handle_multiple_frames_per_burst"),
           D(header_extract,make)
        )
        




        
        .def("get_n_dropped_bursts",&header_extract::get_n_dropped_bursts,       
            D(header_extract,get_n_dropped_bursts)
        )


        
        .def("get_input_queue_size",&header_extract::get_input_queue_size,       
            D(header_extract,get_input_queue_size)
        )


        
        .def("debug_id",&header_extract::debug_id,       
            py::arg("id"),
            D(header_extract,debug_id)
        )

        ;




}








