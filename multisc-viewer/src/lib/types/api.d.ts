export interface paths {
    "/health": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Check daemon health */
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
        /** Get loaded datasets */
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
        /** Unload datasets */
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
        /** Load datasets */
        post: operations["/load-POST"];
        delete?: never;
        options?: never;
        head?: never;
        patch?: never;
        trace?: never;
    };
    "/datasets-index": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        /** Get datasets index */
        get: operations["/datasets-index-GET"];
        put?: never;
        post?: never;
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
        /** Get datasets metadata */
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
        /** Get publications metadata */
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
        /** Get genes for datasets */
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
        /** Get differentially expressed genes for datasets */
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
        /** Query differentially expressed genes across datasets */
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
        /** Generate plots for datasets */
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
            /**
             * @description Dataset id
             * @example example-ds-1
             */
            id: string;
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
            publicationDate?: string;
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
            /**
             * @description Number of genes in dataset
             * @example 20000
             */
            genes?: number;
            deg?: {
                /** @example example-ds-deg */
                id: string;
                /** @example Dataset DEGs */
                name?: string;
                /**
                 * @description Number of differentially expressed genes in dataset
                 * @example 2000
                 */
                genes?: number;
            }[];
        };
        Publication: {
            /** @example example-pub-1 */
            id: string;
            author?: components["schemas"]["Author"];
            /**
             * Format: date
             * @description YYYY-MM-DD
             * @example 2000-01-01
             */
            publicationDate?: string;
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
        Gene: {
            /** @example example-ds-1 */
            id: string;
            /**
             * @description DEG set name
             * @example Dataset DEGs
             */
            name?: string;
            /**
             * @example [
             *       "APOE",
             *       "TREM2",
             *       "CLU"
             *     ]
             */
            gene: string[];
        };
        DEGs: {
            /** @example example-ds-1 */
            id: string;
            deg: components["schemas"]["Gene"];
        };
        GenesRows: {
            /** @example APOE */
            id: string;
            /**
             * @example [
             *       "example-ds-1",
             *       "example-ds-2"
             *     ]
             */
            gene: string[];
            /**
             * @example [
             *       "example-ds-deg1",
             *       "example-ds-deg2"
             *     ]
             */
            deg: string[];
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
        /** @description Generated plot ids */
        Plot: string[];
        /** @description Error response */
        Error: {
            /**
             * Format: uri
             * @example https://datatracker.ietf.org/doc/html/rfc9110#section-15.5.1
             */
            type?: string;
            /** @example Bad Request */
            title?: string;
            /**
             * Format: int32
             * @example 400
             */
            status: number;
            /** @example Detailed error message. */
            detail?: string;
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
            /** @description OK */
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
            /** @description Error response */
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
            /** @description Loaded dataset ids */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": string[];
                };
            };
            /** @description Error response */
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
                     * @description Dataset ids. Unset to unload all datasets
                     * @example [
                     *       "example-ds-1"
                     *     ]
                     */
                    ds?: string[];
                };
            };
        };
        responses: {
            /** @description Loaded dataset ids */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": string[];
                };
            };
            /** @description Error response */
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
                     * @description Dataset ids
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
            /** @description Loaded dataset ids */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": string[];
                };
            };
            /** @description Error response */
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
    "/datasets-index-GET": {
        parameters: {
            query?: never;
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description List of datasets ids and sizes */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": {
                        /** @example example-ds-1 */
                        id: string;
                        /**
                         * @description Size of dataset in bytes
                         * @example 123456789
                         */
                        size: number;
                    }[];
                };
            };
            /** @description Error response */
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
            /** @description List of datasets metadata */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Dataset"][];
                };
            };
            /** @description Error response */
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
            /** @description List of publications metadata */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Publication"][];
                };
            };
            /** @description Error response */
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
            /** @description List of datasets genes */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Gene"][];
                };
            };
            /** @description Error response */
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
            /** @description List of datasets differentially expressed genes */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["DEGs"][];
                };
            };
            /** @description Error response */
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
                /** @description Genes query filter */
                gene?: string;
                /** @description Results limit */
                limit?: number;
                /** @description Results offset */
                offset?: number;
            };
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description List of genes rows */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["GenesRows"][];
                };
            };
            /** @description Error response */
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
            query: {
                /** @description Dataset ids */
                ds: string[];
                /** @description Gene symbols */
                gene: string[];
                /** @description Plot types */
                pt: ("umap" | "vln" | "feat")[];
                /** @description Group by parameter */
                groupBy: "CellType" | "Genotype";
                /** @description Split by parameter */
                splitBy: "CellType" | "Genotype";
            };
            header?: never;
            path?: never;
            cookie?: never;
        };
        requestBody?: never;
        responses: {
            /** @description Generated plot ids */
            200: {
                headers: {
                    [name: string]: unknown;
                };
                content: {
                    "application/json": components["schemas"]["Plot"];
                };
            };
            /** @description Error response */
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
