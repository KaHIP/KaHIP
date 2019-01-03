/******************************************************************************
 * parse_spac_parameters.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2018
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef KAHIP_PARSE_SPAC_PARAMETERS_H
#define KAHIP_PARSE_SPAC_PARAMETERS_H

#include <string>
#include <sstream>
#include "configuration.h"

struct SpacConfig {
    EdgeWeight infinity;
};

int parse_spac_parameters(int argn, char **argv, PartitionConfig &partition_config, SpacConfig &spac_config,
                          std::string &graph_filename) {

    const char *progname = argv[0];

    struct arg_lit *help = arg_lit0(NULL, "help", "Print help.");
    struct arg_str *filename = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file to partition.");
    struct arg_int *k = arg_int1(NULL, "k", NULL, "Number of blocks to partition the graph.");
    struct arg_int *seed = arg_int0(NULL, "seed", NULL, "Seed to use for PRNG.");
    struct arg_rex *preconfiguration = arg_rex1(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: eco) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
    struct arg_int *infinity = arg_int0(NULL, "infinity", NULL, "Infinity edge weight. Default: 1000");
    struct arg_end *end = arg_end(100);

    void *argtable[] = {
            help, filename, k, seed, preconfiguration, infinity, end
    };

    // Parse arguments.
    int nerrors = arg_parse(argn, argv, argtable);

    // Catch case that help was requested.
    if (help->count > 0) {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        arg_print_glossary(stdout, argtable,"  %-40s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }

    if (nerrors > 0) {
        arg_print_errors(stderr, end, progname);
        printf("Try '%s --help' for more information.\n",progname);
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return 1;
    }

    configuration cfg;
    cfg.standard(partition_config);

    if (k->count > 0) {
        partition_config.k = k->ival[0];
    }

    if (filename->count > 0) {
        graph_filename = filename->sval[0];
    }

    if(preconfiguration->count > 0) {
        if(strcmp("strong", preconfiguration->sval[0]) == 0) {
            cfg.strong(partition_config);
        } else if (strcmp("eco", preconfiguration->sval[0]) == 0) {
            cfg.eco(partition_config);
        } else if (strcmp("fast", preconfiguration->sval[0]) == 0) {
            cfg.fast(partition_config);
        } else if (strcmp("fastsocial", preconfiguration->sval[0]) == 0) {
            cfg.fastsocial(partition_config);
        } else if (strcmp("ecosocial", preconfiguration->sval[0]) == 0) {
            cfg.ecosocial(partition_config);
        } else if (strcmp("strongsocial", preconfiguration->sval[0]) == 0) {
            cfg.strongsocial(partition_config);
        } else {
            fprintf(stderr, "Invalid preconfiguration variant: \"%s\"\n", preconfiguration->sval[0]);
            exit(0);
        }
    }

    if (seed->count > 0) {
        partition_config.seed = seed->ival[0];
    }

    if (infinity->count > 0) {
        spac_config.infinity = infinity->ival[0];
    } else {
        spac_config.infinity = 1000;
    }

    return 0;
}

#endif // KAHIP_PARSE_SPAC_PARAMETERS_H
