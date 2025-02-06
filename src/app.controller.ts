import {
  BadRequestException,
  Controller,
  Get,
  Post,
  Query,
} from "@nestjs/common";
import { AppService, ChartResult, Dataset } from "./app.service";
import { Page } from "./utils/decorators/page.decorator";

@Controller()
export class AppController {
  constructor(private readonly service: AppService) {}

  @Page()
  @Get("/")
  public index(): PageProps {
    return {};
  }

  @Page()
  @Get("/compare")
  public comparison(
    @Query("datasets") datasets: string = "",
    @Query("gene") gene: string = "",
    @Query("groupBy") groupBy: string = "",
    @Query("splitBy") splitBy: string = "",
    @Query("token") token: string | null = null,
  ): CompareProps {
    const knownDatasets = Object.fromEntries(
      this.getDatasets().map((ds) => [ds.name, ds]),
    );

    if (!datasets?.split(",").every((dataset) => dataset in knownDatasets))
      throw new BadRequestException("Unknown dataset requested");

    const genes = this.getGenes(datasets);

    if (!gene) {
      const defaultGene = genes[0];
      if (!defaultGene)
        throw new BadRequestException("No common genes between datasets");
      gene = defaultGene;
    }
    if (!genes.includes(gene))
      throw new BadRequestException("Selected gene is not in all datasets");

    if (!groupBy && !splitBy) {
      groupBy = "Genotype";
      splitBy = "CellType";
    } else if (!groupBy)
      groupBy = splitBy === "Genotype" ? "CellType" : "Genotype";
    else if (!splitBy)
      splitBy = groupBy === "Genotype" ? "CellType" : "Genotype";

    this.service.preload(token, datasets.split(","));

    return {
      order: datasets
        .split(",")
        .sort((a, b) => knownDatasets[a].size - knownDatasets[b].size),
      genes,
      gene,
    };
  }

  @Get("/datasets")
  public getDatasets(): Dataset[] {
    return this.service.getDatasets();
  }

  @Get("/genes")
  public getGenes(@Query("datasets") datasets: string): string[] {
    return this.service.getGenes(datasets.split(","));
  }

  @Get("/plot")
  public async getPlot(
    @Query("dataset") dataset: string,
    @Query("gene") gene: string,
    @Query("groupBy") groupBy: string,
    @Query("splitBy") splitBy: string,
    @Query("token") token: string | null = null,
  ): Promise<ChartResult> {
    return this.service.render(dataset, gene, groupBy, splitBy, token);
  }

  @Post("/preload")
  public async preload(
    @Query("token") token: string | null = null,
    @Query("datasets") datasets: string,
  ): Promise<void> {
    return this.service.preload(token, datasets.split(","));
  }
}
