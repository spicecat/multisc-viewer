import { Controller, Get, Query } from '@nestjs/common';
import { AppService, ChartResult, Dataset } from './app.service';
import { Page } from './utils/decorators/page.decorator';

@Controller()
export class AppController {
	constructor(private readonly service: AppService) {}

	@Page()
	@Get('/')
	public index(): PageProps {
		return {};
	}

	@Get('/datasets')
	public getDatasets(): Dataset[] {
		return this.service.getDatasets();
	}

	@Get('/genes')
	public async getGenes(@Query('datasets') datasets: string): Promise<string[]> {
		return this.service.getGenes(datasets.split(','));
	}

	@Get('/plots')
	public async getPlots(
		@Query('datasets') datasets: string,
		@Query('gene') gene: string,
		@Query('groupBy') groupBy: string,
		@Query('splitBy') splitBy: string
	): Promise<Record<string, ChartResult>> {
		return this.service.generate(datasets.split(','), gene, groupBy, splitBy);
	}
}

