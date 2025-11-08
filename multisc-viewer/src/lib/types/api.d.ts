export interface paths {
    "/health": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Check daemon health. */
        get: operations["/health-GET"];
        put?: never;
        post?: never;
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/loaded": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Get loaded datasets on daemon. */
        get: operations["/loaded-GET"];
        put?: never;
        post?: never;
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/unload": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        get?: never;
        put?: never;
        /** Unload datasets on daemon. */
        post: operations["/unload-POST"];
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/load": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        get?: never;
        put?: never;
        /** Load datasets on daemon. */
        post: operations["/load-POST"];
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/datasets": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Get datasets metadata. */
        get: operations["/datasets-GET"];
        put?: never;
        post?: never;
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/publications": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Get publications metadata. */
        get: operations["/publications-GET"];
        put?: never;
        post?: never;
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/genes": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Get genes for datasets. */
        get: operations["/genes-GET"];
        put?: never;
        post?: never;
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/degs": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Get differentially expressed genes for datasets. */
        get: operations["/degs-GET"];
        put?: never;
        post?: never;
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/genes-rows": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Get differentially expressed genes for datasets. */
        get: operations["/genes-rows-GET"];
        put?: never;
        post?: never;
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/plots": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        get?: never;
        put?: never;
        /** Generate plots for datasets. */
        post: operations["/plots-POST"];
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
}
export type webhooks = Record<string, never>;
export interface components {
    schemas: {
        Dataset: {
            _id: string;
            author?: components["schemas"]["Author"];
            citation?: {
                /** @example PMID:00000000 */
                identifier?: string;
                /** @example Primary Study Title */
                name?: string;
                /**
                 * Format: uri
                 * @example https://pubmed.ncbi.nlm.nih.gov/00000000/
                 */
                url?: string;
            }[];
            /**
             * Format: date
             * @description Formatted as YYYY-MM-DD
             * @example 2000-01-01
             */
            date?: string;
            /**
             * @description Default gene to plot
             * @example APOE
             */
            defaultGene?: string;
            /** @example This is a summary of this dataset. */
            description?: string;
            /** @example Dataset 1 */
            displayName?: string;
            healthCondition?: components["schemas"]["OntologyTerm"];
            /** @example GSE000000 */
            identifier?: string;
            /** @example Dataset Title */
            name?: string;
            species?: components["schemas"]["OntologyTerm"];
            cellType?: components["schemas"]["OntologyTerm"];
            tissue?: components["schemas"]["OntologyTerm"];
            /**
             * Format: uri
             * @example https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE000000
             */
            url?: string;
            /** @description Number of genes */
            gene?: number;
            deg?: {
                [key: string]: {
                    /** @example example-ds-deg */
                    _id: string;
                    /** @description Number of differentially expressed genes */
                    gene: number;
                    /** @example Dataset DEGs */
                    name?: string;
                };
            };
            /** @description Size of dataset in bytes */
            size?: number;
        };
        Publication: {
            /** @example example-pub-1 */
            _id: string;
            author?: components["schemas"]["Author"];
            /**
             * Format: date
             * @description YYYY-MM-DD
             * @example 2000-01-01
             */
            date?: string;
            /** @example This is an abstract of this publication. */
            description?: string;
            /** @example Nature */
            journalName?: string;
            /** @example PMID:00000000 */
            identifier?: string;
            /** @example Primary Study Title */
            name?: string;
            /** @example https://pubmed.ncbi.nlm.nih.gov/00000000/ */
            url?: string;
            /**
             * @example [
             *       "example-ds-1",
             *       "example-ds-2"
             *     ]
             */
            datasets: string[];
        };
        Gene: string[];
        DEGs: {
            [key: string]: {
                /** @example example-ds-deg */
                _id: string;
                gene: components["schemas"]["Gene"];
                /** @example Dataset DEGs */
                name?: string;
            };
        };
        Author: {
            /**
             * @description Surname Initial
             * @example Doe J
             */
            name: string;
        }[];
        OntologyTerm: {
            /** @example AD */
            displayName?: string;
            /** @example mesh:D000544 */
            identifier?: string;
            /** @example Alzheimer Disease */
            name: string;
            /**
             * Format: uri
             * @example http://id.nlm.nih.gov/mesh/D000544
             */
            url?: string;
        }[];
        /** @description List of paths to directories with plots */
        Plot: string[];
        /** @description Error response */
        Error: {
            /** @example 500 - Internal server error */
            error: string;
        };
    };
    responses: never;
    parameters: never;
    requestBodies: never;
    headers: never;
    pathItems: never;
}
export type $defs = Record<string, never>;
export interface operations {
    "/health-GET": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description Daemon health */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": {
                        /** @example ok */
                        status: string;
                    };
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/loaded-GET": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description List of loaded dataset ids */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": string[];
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/unload-POST": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody: {
            content: {
                "application/json": {
                    /**
                     * @description Dataset ids
                     * @example [
                     *       "example-ds-1"
                     *     ]
                     */
                    ds?: string[];
                };
            };
        };
        responses: {
            /** @description List of loaded dataset ids */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": string[];
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/load-POST": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: {
            content: {
                "application/json": {
                    /**
                     * @description Dataset ids. Unset to unload all datasets.
                     * @example [
                     *       "example-ds-1",
                     *       "example-ds-2"
                     *     ]
                     */
                    ds?: string[];
                };
            };
        };
        responses: {
            /** @description List of loaded dataset ids */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": string[];
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/datasets-GET": {
        parameters: {
            query?: {
                /** @description Dataset ids. Unset to get all datasets. */
                ds?: string[];
            };
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description Dictionary mapping dataset id to metadata */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": {
                        [key: string]: components["schemas"]["Dataset"];
                    };
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/publications-GET": {
        parameters: {
            query?: {
                /** @description Publication ids. Unset to get all publications. */
                pub?: string[];
            };
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description Dictionary mapping publication id to metadata */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": {
                        [key: string]: components["schemas"]["Publication"];
                    };
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/genes-GET": {
        parameters: {
            query: {
                /** @description Dataset ids */
                ds: string[];
            };
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description Dictionary mapping dataset id to genes */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": {
                        [key: string]: components["schemas"]["Gene"];
                    };
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/degs-GET": {
        parameters: {
            query: {
                /** @description Dataset ids */
                ds: string[];
            };
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description Dictionary mapping dataset id to dictionary mapping DEG id to differentially expressed genes */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": {
                        [key: string]: components["schemas"]["DEGs"];
                    };
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/genes-rows-GET": {
        parameters: {
            query: {
                /** @description Dataset ids */
                ds: string[];
            };
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description List of gene rows */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": {
                        /** @example APOE */
                        _id: string;
                        /**
                         * @example [
                         *       "example-ds-1",
                         *       "example-ds-2"
                         *     ]
                         */
                        datasets: string[];
                        /**
                         * @example [
                         *       "example-ds-deg"
                         *     ]
                         */
                        degs: string[];
                    }[];
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
    "/plots-POST": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody: {
            content: {
                "application/json": {
                    /**
                     * @description Dataset ids
                     * @example [
                     *       "example-ds-1",
                     *       "example-ds-2"
                     *     ]
                     */
                    ds: string[];
                    gene: components["schemas"]["Gene"];
                    /**
                     * @description Plot types
                     * @default [
                     *       "umap",
                     *       "vln",
                     *       "feat"
                     *     ]
                     * @example [
                     *       "umap",
                     *       "vln",
                     *       "feat"
                     *     ]
                     */
                    pt: ("umap" | "vln" | "feat")[];
                    /**
                     * @default CellType
                     * @example CellType
                     * @enum {string}
                     */
                    groupBy: "CellType" | "Genotype";
                    /**
                     * @default Genotype
                     * @example Genotype
                     * @enum {string}
                     */
                    splitBy: "CellType" | "Genotype";
                };
            };
        };
        responses: {
            /** @description List of paths to directories with plots */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Plot"];
                };
            };
            /** @description Unexpected error */
            default: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Error"];
                };
            };
        };
    };
}
