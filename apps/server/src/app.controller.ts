import {
  BadRequestException,
  Controller,
  Get,
  Body,
  Post,
  Query,
} from "@nestjs/common";
import { DataService } from "./data.service";
import { PlotService } from "./plot.service";
import { DaemonService } from "./daemon.service";
// import { Page } from "./utils/decorators/page.decorator";

@Controller()
export class AppController {
  constructor(
    private readonly dataService: DataService,
    private readonly plotService: PlotService,
    private readonly daemonService: DaemonService
  ) {}

  @Get("/publications")
  getPublications(): Publication[] {
    return [...this.dataService.publications.values()];
  }

  @Get("/datasets")
  getDatasets(): Dataset[] {
    return [...this.dataService.datasets.values()];
  }

  @Get("/genes")
  getGenes(@Query("datasets") datasets: string = ""): string[] {
    return this.dataService.getGenes(datasets.split(","));
  }

  @Post("/load")
  async load(@Body("datasets") datasets: string = ""): Promise<void> {
    return this.daemonService.load(datasets.split(","));
  }

  @Get("/plot")
  async getPlot(
    @Query("dataset") dataset: string,
    @Query("gene") gene: string,
    @Query("groupBy") groupBy: string,
    @Query("splitBy") splitBy: string,
  ): Promise<ChartResult> {
    return this.plotService.render(dataset, gene, groupBy, splitBy);
  }

  // @Page()
  // @Get("/")
  // index(): IndexProps {
  //   const datasets = this.service.getDatasets();
  //   return { datasets, token: randomBytes(32).toString("hex") };
  // }

  // @Page()
  // @Get("/publication")
  // publication() {
  //   return {
  //     publications: this.service.getPublications(),
  //   };
  // }

  // @Page()
  // @Get("/publication/:publicationId")
  // publicationId(@Param("publicationId") publicationId: string) {
  //   const publication = this.service
  //     .getPublications()
  //     .find((publication) => publication.publicationId === publicationId);
  //   if (!publication)
  //     throw new BadRequestException("Unknown publication requested");
  //   return { publication };
  // }

  // @Page()
  // @Get("/compare")
  // comparison(
  //   @Query("datasets") datasets: string = "",
  //   @Query("gene") gene: string = "",
  //   @Query("groupBy") groupBy: string = "",
  //   @Query("splitBy") splitBy: string = "",
  //   @Query("token") token: string | null = null,
  // ): CompareProps {
  //   const knownDatasets = Object.fromEntries(
  //     this.getDatasets().map((ds) => [ds.name, ds]),
  //   );

  //   if (!datasets?.split(",").every((ds) => ds in knownDatasets))
  //     throw new BadRequestException("Unknown dataset requested");

  //   const genes = this.getGenes(datasets);

  //   if (!gene) {
  //     const defaultGene = genes[0];
  //     if (!defaultGene)
  //       throw new BadRequestException("No common genes between datasets");
  //     gene = defaultGene;
  //   }
  //   if (!genes.includes(gene))
  //     throw new BadRequestException("Selected gene is not in all datasets");

  //   if (!groupBy && !splitBy) {
  //     groupBy = "Genotype";
  //     splitBy = "CellType";
  //   } else if (!groupBy)
  //     groupBy = splitBy === "Genotype" ? "CellType" : "Genotype";
  //   else if (!splitBy)
  //     splitBy = groupBy === "Genotype" ? "CellType" : "Genotype";

  //   this.service.preload(token, datasets.split(","));

  //   return {
  //     order: datasets
  //       .split(",")
  //       .sort((a, b) => knownDatasets[a].size - knownDatasets[b].size),
  //     genes,
  //     gene,
  //   };
  // }
}
